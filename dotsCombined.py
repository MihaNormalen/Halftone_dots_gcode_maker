import math
import argparse
import xml.etree.ElementTree as ET
import sys
import random

# Poskus uvoza osnovnega generatorja
try:
    from dotPainter import GCodeBaseGenerator
except ImportError:
    class GCodeBaseGenerator:
        def __init__(self, feed_rate, x_offset, y_offset,
                     dip_location_raw, dip_duration_s,
                     dip_wipe_radius, z_wipe_travel_raw,
                     dip_entry_radius, remove_drops_enabled,
                     z_global_offset_val, z_safe_raw, z_safe_dip_raw):
            self.feed_rate = feed_rate
            self.x_offset = x_offset
            self.y_offset = y_offset
            self.gcode = []
            self.z_global_offset = z_global_offset_val
            self.z_safe = z_safe_raw + self.z_global_offset
            self.z_safe_dip = z_safe_dip_raw + self.z_global_offset
            self.z_wipe_travel = z_wipe_travel_raw + self.z_global_offset
            self.dip_location = (dip_location_raw[0], dip_location_raw[1], dip_location_raw[2] + self.z_global_offset)
            self.dip_duration_s = dip_duration_s
            self.dip_wipe_radius = dip_wipe_radius
            self.dip_entry_radius = dip_entry_radius
            self.remove_drops_enabled = remove_drops_enabled
            self._initial_dip_performed = False

        def remove_drops(self, tray_x, tray_y, x, y):
            dx, dy = x - tray_x, y - tray_y
            base_angle = math.atan2(dy, dx)
            angle = base_angle + random.uniform(-0.2, 0.2)
            ux, uy = math.cos(angle), math.sin(angle)
            x1, y1 = tray_x + ux * self.dip_entry_radius, tray_y + uy * self.dip_entry_radius
            x2, y2 = tray_x + ux * self.dip_wipe_radius, tray_y + uy * self.dip_wipe_radius
            self.gcode.append(f"G0 X{x1:.3f} Y{y1:.3f}\nG0 Z{self.z_wipe_travel:.3f}\nG1 X{x2:.3f} Y{y2:.3f} F600")

        def _perform_dip(self, target_x=None, target_y=None):
            self.gcode.append(f"\n; --- Namakanje (Spiralno) ---")
            self.gcode.append(f"G0 Z{self.z_safe_dip:.3f}\nG0 X{self.dip_location[0]:.3f} Y{self.dip_location[1]:.3f}")
            self.gcode.append(f"G1 Z{self.dip_location[2]:.3f} F800")
            d_theta = 0.2
            for i in range(32):
                theta = i * d_theta
                r = (theta / (2*math.pi)) * self.dip_entry_radius
                self.gcode.append(f"G1 X{self.dip_location[0] + r*math.cos(theta):.3f} Y{self.dip_location[1] + r*math.sin(theta):.3f} F1000")
            if self.dip_duration_s > 0: self.gcode.append(f"G4 P{int(self.dip_duration_s * 1000)}")
            self.gcode.append(f"G1 Z{self.dip_location[2]+2:.3f} F800")
            if self.remove_drops_enabled and target_x: self.remove_drops(self.dip_location[0], self.dip_location[1], target_x, target_y)
            self.gcode.append(f"G1 Z{self.z_safe_dip:.3f} F600")

        def save(self, filename):
            with open(filename, 'w') as f: f.write("\n".join(self.gcode))

class SVGCircleGenerator(GCodeBaseGenerator):
    def __init__(self, brush_width, overlap, max_distance_per_dip, fill_type, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.brush_width = brush_width
        self.overlap = overlap
        self.max_distance_per_dip = max_distance_per_dip
        self.fill_type = fill_type 
        self.distance_since_dip = 0

    def calculate_path_length(self, r):
        step = self.brush_width * (1.0 - self.overlap)
        if self.fill_type == 'spiral':
            return (math.pi * (r**2)) / (self.brush_width * (1.0 - self.overlap)) if r > 0 else 0
        else:
            total_len, curr_r = 0, r
            while curr_r > 0:
                total_len += 2 * math.pi * curr_r
                curr_r -= step
            return max(total_len, 0.1)

    def paint_spiral(self, cx, cy, r):
        
        step_distance = self.brush_width * (1.0 - self.overlap)
        theta = 0
        self.gcode.append(f"G0 X{cx:.3f} Y{cy:.3f} Z{self.z_safe:.3f}")
        self.gcode.append(f"G1 Z{self.z_global_offset:.3f} F600")
        while True:
            curr_r = (step_distance * theta) / (2 * math.pi)
            if curr_r > r: break
            self.gcode.append(f"G1 X{cx + curr_r*math.cos(theta):.3f} Y{cy + curr_r*math.sin(theta):.3f} F{self.feed_rate}")
            theta += 0.2
        self.gcode.append(f"G1 Z{self.z_safe:.3f} F800")

    def paint_concentric(self, cx, cy, r):
        
        step = self.brush_width * (1.0 - self.overlap)
        curr_r = r
        first = True
        while curr_r >= 0:
            sx, sy = cx + curr_r, cy
            if first:
                self.gcode.append(f"G0 X{sx:.3f} Y{sy:.3f} Z{self.z_safe:.3f}")
                self.gcode.append(f"G1 Z{self.z_global_offset:.3f} F600")
                first = False
            else:
                self.gcode.append(f"G1 X{sx:.3f} Y{sy:.3f} F{self.feed_rate}")
            if curr_r > 0.05:
                self.gcode.append(f"G2 X{sx:.3f} Y{sy:.3f} I{-curr_r:.3f} J0 F{self.feed_rate}")
            else:
                self.gcode.append(f"G1 X{cx:.3f} Y{cy:.3f} F{self.feed_rate}")
                break
            curr_r -= step
        self.gcode.append(f"G1 Z{self.z_safe:.3f} F800")

    def parse_and_scale_svg(self, svg_path, target_width_mm):
        tree = ET.parse(svg_path)
        root = tree.getroot()
        circles = []
        for c in root.iter():
            if 'circle' in c.tag:
                circles.append({'cx': float(c.get('cx',0)), 'cy': float(c.get('cy',0)), 'r': float(c.get('r',0))})
        if not circles: return []
        min_x = min(c['cx'] - c['r'] for c in circles)
        max_x = max(c['cx'] + c['r'] for c in circles)
        max_y = max(c['cy'] + c['r'] for c in circles)
        scale = target_width_mm / (max_x - min_x) if (max_x - min_x) > 0 else 1.0
        return [{'x': (c['cx'] - min_x) * scale + self.x_offset, 'y': (max_y - c['cy']) * scale + self.y_offset, 'r': c['r'] * scale} for c in circles]

    def run(self, svg_path, target_width):
        self.gcode.append("G90\nG21")
        circles = self.parse_and_scale_svg(svg_path, target_width)
        current_pos = (self.dip_location[0], self.dip_location[1])
        
        while circles:
            best_idx = min(range(len(circles)), key=lambda i: math.hypot(circles[i]['x']-current_pos[0], circles[i]['y']-current_pos[1]))
            dot = circles.pop(best_idx)
            path_len = self.calculate_path_length(dot['r'])
            
            if (self.distance_since_dip + path_len) > self.max_distance_per_dip:
                self._perform_dip(target_x=dot['x'], target_y=dot['y'])
                self.distance_since_dip = 0

            if self.fill_type == 'spiral':
                self.paint_spiral(dot['x'], dot['y'], dot['r'])
            else:
                self.paint_concentric(dot['x'], dot['y'], dot['r'])
            
            self.distance_since_dip += path_len
            current_pos = (dot['x'], dot['y'])
        self.gcode.append("M2")

if __name__ == "__main__":
    # ==========================================
    # GLAVNE NASTAVITVE (STOLPEC SPREMENLJIVK)
    # ==========================================
    FILL_TYPE    = 'spiral'    # Izbira: 'spiral' ali 'concentric'
    BRUSH_WIDTH  = 1.5         # Debelina čopiča v mm
    OVERLAP      = 0.8         # Prekrivanje (0.0 - 1.0)
    MAX_DIST     = 1000.0      # Max pot barvanja pred namakanjem (mm)
    TARGET_WIDTH = 200.0       # Ciljna širina slike (mm)
    FEED_RATE    = 1500        # Hitrost premikanja (mm/min)
    
    # Lokacija namakanja in offseti
    X_OFFSET     = 33.0
    Y_OFFSET     = 110.0
    DIP_X        = 91.0
    DIP_Y        = 25.0
    DIP_Z        = 1.3
    DIP_WIPE_R   = 29.0
    
    # Višine Z
    Z_SAFE       = 5.0
    Z_SAFE_DIP   = 16.0
    Z_GLOBAL_OFF = 0.0
    # ==========================================

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Vhodna SVG datoteka")
    parser.add_argument("output", help="Izhodna G-koda")
    args = parser.parse_args()

    painter = SVGCircleGenerator(
        brush_width=BRUSH_WIDTH, overlap=OVERLAP, max_distance_per_dip=MAX_DIST,
        fill_type=FILL_TYPE, feed_rate=FEED_RATE, x_offset=X_OFFSET, y_offset=Y_OFFSET,
        dip_location_raw=(DIP_X, DIP_Y, DIP_Z), dip_duration_s=0.2,
        dip_wipe_radius=DIP_WIPE_R, z_wipe_travel_raw=1.0, dip_entry_radius=3.0,
        remove_drops_enabled=True, z_global_offset_val=Z_GLOBAL_OFF,
        z_safe_raw=Z_SAFE, z_safe_dip_raw=Z_SAFE_DIP
    )
    painter.run(args.input, TARGET_WIDTH)
    painter.save(args.output)
