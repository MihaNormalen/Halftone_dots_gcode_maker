import math
import argparse
import xml.etree.ElementTree as ET
import sys
import random


class GCodeBaseGenerator:
    def __init__(self, feed_travel, feed_paint, accel_travel, accel_paint,
                 x_offset, y_offset,
                 dip_location_raw, dip_jitter, dip_spiral_loops, dip_spiral_r,
                 dip_wipe_radius, z_wipe_exit_raw,
                 z_global_offset_val, z_safe_raw, z_safe_dip_raw,
                 min_dist, max_dist):
        self.feed_travel      = feed_travel
        self.feed_paint       = feed_paint
        self.accel_travel     = accel_travel
        self.accel_paint      = accel_paint
        self.x_offset         = x_offset
        self.y_offset         = y_offset
        self.gcode            = []

        self.z_global_offset  = z_global_offset_val
        self.z_safe           = z_safe_raw     + self.z_global_offset   # high travel Z
        self.z_safe_dip       = z_safe_dip_raw + self.z_global_offset   # Z before/after dip
        self.z_low            = z_safe_raw     + self.z_global_offset   # same as z_safe; lifted but low
        self.z_wipe_exit      = z_wipe_exit_raw + self.z_global_offset
        self.z_paint          = self.z_global_offset                     # Z while painting (0 + offset)

        dip_x, dip_y, dip_z  = dip_location_raw
        self.dip_x            = dip_x
        self.dip_y            = dip_y
        self.dip_z            = dip_z + self.z_global_offset
        self.dip_jitter       = dip_jitter
        self.dip_spiral_loops = dip_spiral_loops
        self.dip_spiral_r     = dip_spiral_r
        self.dip_wipe_radius  = dip_wipe_radius

        self.min_dist         = min_dist
        self.max_dist         = max_dist
        self.dist_since_dip   = 0.0
        self.current_max_dist = random.uniform(min_dist, max_dist)
        self.current_pos      = (dip_x, dip_y)

    # ------------------------------------------------------------------
    # Speed helpers
    # ------------------------------------------------------------------

    def _set_speed(self, mode='travel'):
        """Switch acceleration + feed rate.  M400 flushes the planner first."""
        accel = self.accel_travel if mode == 'travel' else self.accel_paint
        feed  = self.feed_travel  if mode == 'travel' else self.feed_paint
        self.gcode.append("M400")
        self.gcode.append(f"M204 P{accel} T{accel}")
        self.gcode.append(f"G1 F{feed}")

    # ------------------------------------------------------------------
    # Dip cycle (ported from lines500)
    # ------------------------------------------------------------------

    def _perform_dip_and_travel(self, target_x, target_y):
        """Full dip cycle: lift → travel to dip (with jitter) → spiral →
        directional wipe → travel to target.  Resets distance counter."""
        self.gcode.append(f"\n; --- CIKEL NAMAKANJA ---")

        # Lift to safe Z before moving
        self.gcode.append(f"G0 Z{self.z_safe_dip:.3f} F3000")
        self._set_speed('travel')

        # Dip position with random jitter
        ax = self.dip_x + random.uniform(-self.dip_jitter, self.dip_jitter)
        ay = self.dip_y + random.uniform(-self.dip_jitter, self.dip_jitter)

        # Fly to dip XY at safe Z, then plunge
        self.gcode.append(f"G0 X{ax:.3f} Y{ay:.3f} Z{self.z_safe_dip:.3f}")
        self.gcode.append(f"G1 Z{self.dip_z:.3f} F3000")

        # Outward spiral while submerged
        num_steps = int(self.dip_spiral_loops * 4)
        for i in range(num_steps):
            ang = i * (math.pi / 2)
            r   = (i / max(num_steps - 1, 1)) * self.dip_spiral_r
            self.gcode.append(
                f"G1 X{ax + r*math.cos(ang):.3f} Y{ay + r*math.sin(ang):.3f} F2500"
            )

        # Directional wipe: wipe in the direction of the next target
        dx, dy = target_x - self.dip_x, target_y - self.dip_y
        dist   = math.hypot(dx, dy)
        if dist > 0:
            wx = self.dip_x + (dx / dist) * self.dip_wipe_radius
            wy = self.dip_y + (dy / dist) * self.dip_wipe_radius
        else:
            wx = self.dip_x + self.dip_wipe_radius
            wy = self.dip_y

        self.gcode.append(f"G0 Z{self.z_wipe_exit:.3f} F3000")
        self.gcode.append(f"G0 X{wx:.3f} Y{wy:.3f} F{self.feed_travel}")
        self.gcode.append(f"G0 Z{self.z_safe_dip:.3f} F3000")

        # Travel to target at safe Z then descend to low Z
        self.gcode.append(f"G0 X{target_x:.3f} Y{target_y:.3f} Z{self.z_low:.3f} F{self.feed_travel}")

        # Reset dip counters
        self.dist_since_dip   = 0.0
        self.current_max_dist = random.uniform(self.min_dist, self.max_dist)
        self.current_pos      = (target_x, target_y)

    # ------------------------------------------------------------------

    def save(self, filename):
        with open(filename, 'w') as f:
            f.write("\n".join(self.gcode))


# ======================================================================

class SVGCircleGenerator(GCodeBaseGenerator):
    def __init__(self, brush_width, overlap, fill_type, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.brush_width = brush_width
        self.overlap     = overlap
        self.fill_type   = fill_type

    # ------------------------------------------------------------------
    # Path-length estimation (for dip budget)
    # ------------------------------------------------------------------

    def calculate_path_length(self, r):
        step = self.brush_width * (1.0 - self.overlap)
        if self.fill_type == 'spiral':
            return (math.pi * r ** 2) / step if r > 0 else 0
        else:
            total, curr_r = 0.0, r
            while curr_r > 0:
                total  += 2 * math.pi * curr_r
                curr_r -= step
            return max(total, 0.1)

    # ------------------------------------------------------------------
    # Paint routines
    # ------------------------------------------------------------------

    def _travel_to(self, x, y):
        """G0 travel move to XY at z_low, with speed set to travel mode."""
        self._set_speed('travel')
        self.gcode.append(f"G0 X{x:.3f} Y{y:.3f} Z{self.z_low:.3f}")

    def _check_and_dip(self, needed_dist, target_x, target_y):
        """Trigger a dip cycle if the dip budget would be exceeded."""
        if (self.dist_since_dip + needed_dist) > self.current_max_dist:
            self.gcode.append(f"G0 Z{self.z_low:.3f} F3000")
            self._perform_dip_and_travel(target_x, target_y)

    def paint_spiral(self, cx, cy, r):
        """Archimedean spiral fill, with mid-paint dip support."""
        step_distance = self.brush_width * (1.0 - self.overlap)
        path_len      = self.calculate_path_length(r)

        # Dip before starting if budget exceeded
        self._check_and_dip(path_len, cx, cy)

        # Travel to circle centre
        self._travel_to(cx, cy)
        self._set_speed('paint')
        self.gcode.append(f"G1 Z{self.z_paint:.3f} F600")

        theta = 0.0
        while True:
            curr_r = (step_distance * theta) / (2 * math.pi)
            if curr_r > r:
                break
            px = cx + curr_r * math.cos(theta)
            py = cy + curr_r * math.sin(theta)

            seg_len = math.hypot(px - self.current_pos[0], py - self.current_pos[1])
            # Mid-spiral dip if needed
            if (self.dist_since_dip + seg_len) > self.current_max_dist:
                self.gcode.append(f"G0 Z{self.z_low:.3f} F3000")
                self._perform_dip_and_travel(px, py)
                self._set_speed('paint')
                self.gcode.append(f"G1 Z{self.z_paint:.3f} F600")

            self.gcode.append(f"G1 X{px:.3f} Y{py:.3f} F{self.feed_paint}")
            self.dist_since_dip += seg_len
            self.current_pos     = (px, py)
            theta += 0.2

        # Lift after dot
        self.gcode.append(f"G0 Z{self.z_low:.3f} F3000")
        self.current_pos = (cx, cy)

    def paint_concentric(self, cx, cy, r):
        """Concentric circles fill, outer → inner, with mid-paint dip support."""
        step   = self.brush_width * (1.0 - self.overlap)
        curr_r = r
        first  = True

        while curr_r >= 0:
            sx, sy   = cx + curr_r, cy
            ring_len = 2 * math.pi * curr_r if curr_r > 0.05 else 0

            # Dip check before each ring
            self._check_and_dip(ring_len, sx, sy)

            if first:
                self._travel_to(sx, sy)
                self._set_speed('paint')
                self.gcode.append(f"G1 Z{self.z_paint:.3f} F600")
                first = False
            else:
                # Travel to ring start in paint mode (already at paint Z)
                self._set_speed('travel')
                self.gcode.append(f"G0 X{sx:.3f} Y{sy:.3f}")
                self._set_speed('paint')

            if curr_r > 0.05:
                self.gcode.append(
                    f"G2 X{sx:.3f} Y{sy:.3f} I{-curr_r:.3f} J0 F{self.feed_paint}"
                )
                self.dist_since_dip += ring_len
                self.current_pos     = (sx, sy)
            else:
                self.gcode.append(f"G1 X{cx:.3f} Y{cy:.3f} F{self.feed_paint}")
                self.dist_since_dip += math.hypot(sx - cx, sy - cy)
                self.current_pos     = (cx, cy)
                break

            curr_r -= step

        # Lift after dot
        self.gcode.append(f"G0 Z{self.z_low:.3f} F3000")

    # ------------------------------------------------------------------
    # SVG parsing
    # ------------------------------------------------------------------

    def parse_and_scale_svg(self, svg_path, target_width_mm):
        tree    = ET.parse(svg_path)
        root    = tree.getroot()
        circles = []
        for c in root.iter():
            if 'circle' in c.tag:
                circles.append({
                    'cx': float(c.get('cx', 0)),
                    'cy': float(c.get('cy', 0)),
                    'r':  float(c.get('r',  0)),
                })
        if not circles:
            return []
        min_x = min(c['cx'] - c['r'] for c in circles)
        max_x = max(c['cx'] + c['r'] for c in circles)
        max_y = max(c['cy'] + c['r'] for c in circles)
        scale = target_width_mm / (max_x - min_x) if (max_x - min_x) > 0 else 1.0
        return [
            {
                'x': (c['cx'] - min_x) * scale + self.x_offset,
                'y': (max_y - c['cy']) * scale + self.y_offset,
                'r': c['r'] * scale,
            }
            for c in circles
        ]

    # ------------------------------------------------------------------
    # Main run loop — nearest-neighbour traversal
    # ------------------------------------------------------------------

    def run(self, svg_path, target_width):
        self.gcode.append("G90\nG21")
        circles = self.parse_and_scale_svg(svg_path, target_width)

        # Initial dip before any painting
        if circles:
            first = min(circles, key=lambda c: math.hypot(
                c['x'] - self.dip_x, c['y'] - self.dip_y))
            self._perform_dip_and_travel(first['x'], first['y'])

        while circles:
            # Nearest-neighbour pick
            best_idx = min(
                range(len(circles)),
                key=lambda i: math.hypot(
                    circles[i]['x'] - self.current_pos[0],
                    circles[i]['y'] - self.current_pos[1],
                ),
            )
            dot = circles.pop(best_idx)

            # Travel to dot start (G0, travel speed)
            self._travel_to(dot['x'], dot['y'])

            if self.fill_type == 'spiral':
                self.paint_spiral(dot['x'], dot['y'], dot['r'])
            else:
                self.paint_concentric(dot['x'], dot['y'], dot['r'])

            self.current_pos = (dot['x'], dot['y'])

        self.gcode.append("M400")
        self.gcode.append(f"G0 Z{self.z_safe_dip:.3f} F3000\nM2")


# ======================================================================

if __name__ == "__main__":
    # ==========================================
    # GLAVNE NASTAVITVE (STOLPEC SPREMENLJIVK)
    # ==========================================
    FILL_TYPE    = 'spiral'    # 'spiral' ali 'concentric'
    BRUSH_WIDTH  = 1.1         # Debelina čopiča (mm)
    OVERLAP      = 0.1         # Prekrivanje (0.0 – 1.0)
    TARGET_WIDTH = 270.0       # Ciljna širina slike (mm)

    # Hitrosti & pospeški
    FEED_TRAVEL  = 12000       # Potovalna hitrost (mm/min)
    FEED_PAINT   = 600         # Hitrost barvanja (mm/min)
    ACCEL_TRAVEL = 12000       # Pospešek pri potovanju (mm/s²)
    ACCEL_PAINT  = 1000        # Pospešek pri barvanju (mm/s²)

    # Doseg čopiča med dopi
    MIN_DIST     = 140.0       # Min pot pred dip-om (mm)
    MAX_DIST     = 180.0       # Max pot pred dip-om (mm)

    # Lokacija namakanja in offseti
    X_OFFSET     = 10.0
    Y_OFFSET     = 50.0
    DIP_X        = 88.0
    DIP_Y        = 0.0
    DIP_Z        = 1.3         # Z globina namakanja
    DIP_JITTER   = 5.0        # Naključni odmik namakanja (mm)
    DIP_SPIRAL_LOOPS = 1.0     # Stevilo zavojev spirale pri dipu
    DIP_SPIRAL_R = 5.0        # Radij spirale pri dipu (mm)
    DIP_WIPE_R   = 32.2        # Razdalja brisanja po dipu (mm)

    # Višine Z
    Z_SAFE       = 1.6         # Nizka varna Z (nad površino)
    Z_SAFE_DIP   = 8.0        # Visoka varna Z (nad vsem)
    Z_WIPE_EXIT  = 4.4         # Z pri izhodu iz brisanja
    Z_GLOBAL_OFF = 0.0
    # ==========================================

    parser = argparse.ArgumentParser()
    parser.add_argument("input",  help="Vhodna SVG datoteka")
    parser.add_argument("output", help="Izhodna G-koda")
    args = parser.parse_args()

    painter = SVGCircleGenerator(
        brush_width=BRUSH_WIDTH,
        overlap=OVERLAP,
        fill_type=FILL_TYPE,
        # speed
        feed_travel=FEED_TRAVEL,
        feed_paint=FEED_PAINT,
        accel_travel=ACCEL_TRAVEL,
        accel_paint=ACCEL_PAINT,
        # position
        x_offset=X_OFFSET,
        y_offset=Y_OFFSET,
        # dip geometry
        dip_location_raw=(DIP_X, DIP_Y, DIP_Z),
        dip_jitter=DIP_JITTER,
        dip_spiral_loops=DIP_SPIRAL_LOOPS,
        dip_spiral_r=DIP_SPIRAL_R,
        dip_wipe_radius=DIP_WIPE_R,
        z_wipe_exit_raw=Z_WIPE_EXIT,
        # Z levels
        z_global_offset_val=Z_GLOBAL_OFF,
        z_safe_raw=Z_SAFE,
        z_safe_dip_raw=Z_SAFE_DIP,
        # dip distance budget
        min_dist=MIN_DIST,
        max_dist=MAX_DIST,
    )
    painter.run(args.input, TARGET_WIDTH)
    painter.save(args.output)
