import math
import argparse
import xml.etree.ElementTree as ET
import sys
import random

class GCodeBaseGenerator:
    def __init__(self, feed_rate, feed_rate_paint, x_offset, y_offset,
                 dip_location_raw, dip_duration_s,
                 dip_wipe_radius, z_wipe_travel_raw,
                 dip_entry_radius, remove_drops_enabled,
                 z_global_offset_val, z_safe_raw, z_safe_dip_raw,
                 z_wipe_exit_raw, # New parameter for the low-height exit travel
                 accel_travel=200, accel_paint=20):
        self.feed_rate = feed_rate 
        self.feed_rate_paint = feed_rate_paint 
        self.accel_travel = accel_travel
        self.accel_paint = accel_paint
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.gcode = []
        self.z_global_offset = z_global_offset_val
        self.z_safe = z_safe_raw + self.z_global_offset
        self.z_safe_dip = z_safe_dip_raw + self.z_global_offset
        self.z_wipe_exit = z_wipe_exit_raw + self.z_global_offset
        self.z_wipe_travel = z_wipe_travel_raw + self.z_global_offset
        self.dip_location = (dip_location_raw[0], dip_location_raw[1], dip_location_raw[2] + self.z_global_offset)
        self.dip_duration_s = dip_duration_s
        self.dip_wipe_radius = dip_wipe_radius
        self.dip_entry_radius = dip_entry_radius
        self.remove_drops_enabled = remove_drops_enabled
        self.dip_count = 0
        self.current_pos = (0, 0)

    def _set_machine_speed(self, speed_type='travel'):
        if speed_type == 'travel':
            f = self.feed_rate
            accel = self.accel_travel
        else:
            f = self.feed_rate_paint
            accel = self.accel_paint
        self.gcode.append("M400") 
        self.gcode.append(f"M204 P{accel} T{accel}") 
        self.gcode.append(f"M203 X{f} Y{f} Z{f}") 
        self.gcode.append(f"G0 F{f}") 

    def _update_pos(self, x, y):
        self.current_pos = (x, y)

    def _perform_dip(self, target_x=None, target_y=None):
        self.gcode.append(f"\n; --- Namakanje (Wipe-on-Exit Logic) ---")
        
        j_x = random.uniform(-3.0, 3.0)
        j_y = random.uniform(-3.0, 3.0)
        active_x = self.dip_location[0] + j_x
        active_y = self.dip_location[1] + j_y
        
        # --- 1. DIAGONAL IN ---
        self._set_machine_speed('travel')
        dx = self.dip_location[0] - self.current_pos[0]
        dy = self.dip_location[1] - self.current_pos[1]
        dist = math.hypot(dx, dy)
        
        if dist > 0:
            dx, dy = dx / dist, dy / dist
            edge_x_in = self.dip_location[0] - dx * self.dip_wipe_radius
            edge_y_in = self.dip_location[1] - dy * self.dip_wipe_radius
        else:
            edge_x_in, edge_y_in = self.dip_location[0] + self.dip_wipe_radius, self.dip_location[1]
        
        self.gcode.append(f"G0 X{edge_x_in:.3f} Y{edge_y_in:.3f} Z{self.z_safe_dip:.3f}")
        self.gcode.append(f"G0 X{active_x:.3f} Y{active_y:.3f}")
        
        # --- 2. MIXING ---
        self._set_machine_speed('paint')
        self.gcode.append(f"G1 Z{self.dip_location[2]:.3f} F800")
        direction = 1 if (self.dip_count % 2 == 0) else -1
        self.dip_count += 1
        theta = 0
        max_theta = 2.5 * math.pi
        while theta <= max_theta:
            r = (theta / max_theta) * self.dip_entry_radius
            self.gcode.append(f"G1 X{active_x + r*math.cos(theta*direction):.3f} Y{active_y + r*math.sin(theta*direction):.3f} F1000")
            theta += 0.1
        if self.dip_duration_s > 0: 
            self.gcode.append(f"G4 P{int(self.dip_duration_s * 1000)}")

        # --- 3. WIPE-ON-EXIT LOGIC ---
        self._set_machine_speed('travel')
        
        if target_x is not None and target_y is not None:
            # Calculate exit edge point toward the target
            dx_out = target_x - self.dip_location[0]
            dy_out = target_y - self.dip_location[1]
            dist_out = math.hypot(dx_out, dy_out)
            
            if dist_out > 0:
                dx_out, dy_out = dx_out / dist_out, dy_out / dist_out
                wipe_edge_x = self.dip_location[0] + dx_out * self.dip_wipe_radius
                wipe_edge_y = self.dip_location[1] + dy_out * self.dip_wipe_radius
            else:
                wipe_edge_x, wipe_edge_y = self.dip_location[0] + self.dip_wipe_radius, self.dip_location[1]
            
            # A: Lift to LOW wipe height while inside
            self.gcode.append(f"G0 Z{self.z_wipe_exit:.3f}")
            
            # B: Travel to EDGE at LOW height (The Wipe)
            self.gcode.append(f"G0 X{wipe_edge_x:.3f} Y{wipe_edge_y:.3f}")
            
            # C: Diagonal traverse OUT to target (Climb to Z_SAFE_DIP then drop to Z_SAFE)
            # We first reach max height to clear the dish rim
            self.gcode.append(f"G0 Z{self.z_safe_dip:.3f}")
            self.gcode.append(f"G0 X{target_x:.3f} Y{target_y:.3f} Z{self.z_safe:.3f}")
            
            self._update_pos(target_x, target_y)
        else:
            self.gcode.append(f"G0 Z{self.z_safe_dip:.3f}")

    def save(self, filename):
        with open(filename, 'w') as f: f.write("\n".join(self.gcode))

class SVGCircleGenerator(GCodeBaseGenerator):
    def __init__(self, brush_width, overlap, min_distance_per_dip, max_distance_per_dip, fill_type, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.brush_width = brush_width
        self.overlap = overlap
        self.min_distance_per_dip = min_distance_per_dip
        self.max_distance_per_dip = max_distance_per_dip
        self.fill_type = fill_type 
        self.distance_since_dip = 0
        self.current_max_dist = random.uniform(min_distance_per_dip, max_distance_per_dip)
        self.current_pos = (0, 0)

    def calculate_path_length(self, r):
        step = self.brush_width * (1.0 - self.overlap)
        if self.fill_type == 'spiral':
            return (math.pi * (r**2)) / step if r > 0 else 0
        else:
            total_len, curr_r = 0, r
            while curr_r > 0:
                total_len += 2 * math.pi * curr_r
                curr_r -= step
            return max(total_len, 0.1)

    def paint_spiral(self, cx, cy, r):
        step_distance = self.brush_width * (1.0 - self.overlap)
        theta = 0
        self._set_machine_speed('travel')
        self.gcode.append(f"G0 X{cx:.3f} Y{cy:.3f} Z{self.z_safe:.3f}")
        self._set_machine_speed('paint')
        self.gcode.append(f"G1 Z{self.z_global_offset:.3f} F600")
        while True:
            curr_r = (step_distance * theta) / (2 * math.pi)
            if curr_r > r: break
            self.gcode.append(f"G1 X{cx + curr_r*math.cos(theta):.3f} Y{cy + curr_r*math.sin(theta):.3f} F{self.feed_rate_paint}")
            theta += 0.2
        self._set_machine_speed('travel')
        self.gcode.append(f"G0 Z{self.z_safe:.3f}")
        self._update_pos(cx, cy)

    def paint_concentric(self, cx, cy, r):
        step = self.brush_width * (1.0 - self.overlap)
        curr_r = r
        first = True
        self._set_machine_speed('travel')
        while curr_r >= 0:
            sx, sy = cx + curr_r, cy
            if first:
                self.gcode.append(f"G0 X{sx:.3f} Y{sy:.3f} Z{self.z_safe:.3f}")
                self._set_machine_speed('paint')
                self.gcode.append(f"G1 Z{self.z_global_offset:.3f} F600")
                first = False
            else:
                self.gcode.append(f"G1 X{sx:.3f} Y{sy:.3f} F{self.feed_rate_paint}")
            if curr_r > 0.05:
                mx, my = cx - curr_r, cy
                self.gcode.append(f"G2 X{mx:.3f} Y{my:.3f} I{-curr_r:.3f} J0 F{self.feed_rate_paint}")
                self.gcode.append(f"G2 X{sx:.3f} Y{sy:.3f} I{curr_r:.3f} J0 F{self.feed_rate_paint}")
            else:
                self.gcode.append(f"G1 X{cx:.3f} Y{cy:.3f} F{self.feed_rate_paint}")
                break
            curr_r -= step
        self._set_machine_speed('travel')
        self.gcode.append(f"G0 Z{self.z_safe:.3f}")
        self._update_pos(cx, cy)

    def parse_and_scale_svg(self, svg_path, target_width_mm):
        tree = ET.parse(svg_path)
        root = tree.getroot()
        circles = []
        for c in root.iter():
            if c.tag.endswith('circle'):
                try:
                    cx = float(c.get('cx', 0)); cy = float(c.get('cy', 0)); r = float(c.get('r', 0))
                    circles.append({'cx': cx, 'cy': cy, 'r': r})
                except ValueError: continue
        if not circles: return []
        min_x = min(c['cx'] - c['r'] for c in circles)
        max_x = max(c['cx'] + c['r'] for c in circles)
        max_y = max(c['cy'] + c['r'] for c in circles)
        width = max_x - min_x
        scale = target_width_mm / width if width > 0 else 1.0
        return [{'x': (c['cx'] - min_x) * scale + self.x_offset, 
                 'y': (max_y - c['cy']) * scale + self.y_offset, 
                 'r': c['r'] * scale} for c in circles]

    def run(self, svg_path, target_width):
        self.gcode.append("G90\nG21")
        self._set_machine_speed('travel')
        self.gcode.append(f"G0 Z{self.z_safe_dip:.3f}")
        circles = self.parse_and_scale_svg(svg_path, target_width)
        if not circles: return
        grid_size = 40.0
        grid = {}
        for i, c in enumerate(circles):
            key = (int(c['x'] / grid_size), int(c['y'] / grid_size))
            if key not in grid: grid[key] = []
            grid[key].append(i)
        active = set(range(len(circles)))
        self._perform_dip()
        self.distance_since_dip = 0
        while active:
            best_idx = -1
            min_dist = float('inf')
            gx, gy = int(self.current_pos[0] / grid_size), int(self.current_pos[1] / grid_size)
            found = False
            for radius in range(0, 30):
                if found: break
                for dx in range(-radius, radius + 1):
                    for dy in range(-radius, radius + 1):
                        if abs(dx) != radius and abs(dy) != radius: continue
                        cell = (gx + dx, gy + dy)
                        if cell in grid:
                            for idx in grid[cell]:
                                if idx in active:
                                    d = math.hypot(circles[idx]['x'] - self.current_pos[0], circles[idx]['y'] - self.current_pos[1])
                                    if d < min_dist: min_dist = d; best_idx = idx; found = True
            if best_idx == -1:
                if active: best_idx = min(active, key=lambda i: math.hypot(circles[i]['x'] - self.current_pos[0], circles[i]['y'] - self.current_pos[1]))
                else: break
            active.remove(best_idx)
            dot = circles[best_idx]
            path_len = self.calculate_path_length(dot['r'])
            if (self.distance_since_dip + path_len) > self.current_max_dist:
                self._perform_dip(target_x=dot['x'], target_y=dot['y'])
                self.distance_since_dip = 0
                self.current_max_dist = random.uniform(self.min_distance_per_dip, self.max_distance_per_dip)
            if self.fill_type == 'spiral': self.paint_spiral(dot['x'], dot['y'], dot['r'])
            else: self.paint_concentric(dot['x'], dot['y'], dot['r'])
            self.distance_since_dip += path_len
        self.gcode.append("M2")

if __name__ == "__main__":
    # ==========================================
    # SETTINGS
    # ==========================================
    FILL_TYPE    = 'concentric' 
    BRUSH_WIDTH  = 1.0         
    OVERLAP      = 0.2         
    MIN_DIST     = 1000.0       
    MAX_DIST     = 1500.0       
    TARGET_WIDTH = 300.0      
    FEED_RATE       = 1500     
    FEED_RATE_PAINT = 300      
    ACCEL_TRAVEL    = 200      
    ACCEL_PAINT     = 20       
    X_OFFSET     = 55.0
    Y_OFFSET     = 50.0
    DIP_X        = 15.0
    DIP_Y        = 19.0
    DIP_Z        = 2.3
    DIP_WIPE_R   = 15.0
    
    # HEIGHTS
    Z_SAFE       = 1.8     # Surface travel height
    Z_WIPE_EXIT  = 13.0     # <--- HEIGHT FOR TRAVELING TO RIM (Low exit)
    Z_SAFE_DIP   = 17.0    # High travel height
    Z_GLOBAL_OFF = 0.0
    # ==========================================

    parser = argparse.ArgumentParser()
    parser.add_argument("input"); parser.add_argument("output")
    args = parser.parse_args()

    painter = SVGCircleGenerator(
        brush_width=BRUSH_WIDTH, overlap=OVERLAP, 
        min_distance_per_dip=MIN_DIST, max_distance_per_dip=MAX_DIST,
        fill_type=FILL_TYPE, feed_rate=FEED_RATE, feed_rate_paint=FEED_RATE_PAINT,
        x_offset=X_OFFSET, y_offset=Y_OFFSET,
        dip_location_raw=(DIP_X, DIP_Y, DIP_Z), dip_duration_s=0.2,
        dip_wipe_radius=DIP_WIPE_R, z_wipe_travel_raw=5.0, dip_entry_radius=10.0,
        remove_drops_enabled=True, z_global_offset_val=Z_GLOBAL_OFF,
        z_safe_raw=Z_SAFE, z_safe_dip_raw=Z_SAFE_DIP,
        z_wipe_exit_raw=Z_WIPE_EXIT, # Passed to generator
        accel_travel=ACCEL_TRAVEL, accel_paint=ACCEL_PAINT
    )
    painter.run(args.input, TARGET_WIDTH)
    painter.save(args.output)
