#!/usr/bin/env python3

import argparse
import re
import colorsys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
import numpy as np
from matplotlib.collections import LineCollection


def boost_rgb(rgb_str, saturation_boost=2.0, value_floor=0.85):
    """
    Convert BED RGB string '255,0,0' to a boosted (r, g, b) float tuple.
    - saturation_boost: multiplier on HSV saturation (capped at 1.0)
    - value_floor: minimum brightness — lifts dark colors so they render vividly
    """
    try:
        r, g, b = [int(x) / 255.0 for x in rgb_str.split(',')]
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        s = min(1.0, s * saturation_boost)
        v = max(value_floor, v)
        return colorsys.hsv_to_rgb(h, s, v)
    except Exception:
        return (0.8, 0.8, 0.8)


def boost_rgb_hex(rgb_str, saturation_boost=2.0, value_floor=0.85):
    r, g, b = boost_rgb(rgb_str, saturation_boost, value_floor)
    return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"


class AlignmentVisualizer:
    def __init__(self, bed_files, cigar_files, show_mismatches):
        self.bed_files = bed_files
        self.cigar_files = cigar_files
        self.regions = []
        self.alignments = []
        self.show_mismatches = show_mismatches
        self.sequence_lengths = []
        
    def parse_bed_file(self, bed_file):
        regions = []
        with open(bed_file) as f:
            lines = f.readlines()
            regions = [
                {
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': fields[3],
                    'score': float(fields[4]),
                    'strand': fields[5],
                    'thick_start': int(fields[6]),
                    'thick_end': int(fields[7]),
                    'rgb': fields[8]
                }
                for line in lines
                if len((fields := line.strip().split('\t'))) >= 9
            ]
        return sorted(regions, key=lambda x: x['start'])

    def parse_cigar_string(self, cigar_string):
        pattern = re.compile(r'(\d+)([=XDI])')
        operations = []
        pos1 = pos2 = 0
        
        for length, op in pattern.findall(cigar_string):
            length = int(length)
            operations.append({
                'op': op,
                'length': length,
                'pos1': pos1,
                'pos2': pos2
            })
            
            if op in ['=', 'X']:
                pos1 += length
                pos2 += length
            elif op == 'D':
                pos1 += length
            elif op == 'I':
                pos2 += length
                
        return operations

    def read_cigar_file(self, cigar_file):
        with open(cigar_file) as f:
            return f.readline().strip()

    def process_files(self):
        for bed_file in self.bed_files:
            regions = self.parse_bed_file(bed_file)
            self.regions.append(regions)
            self.sequence_lengths.append(max(region['end'] for region in regions) if regions else 0)
        
        self.alignments = [self.parse_cigar_string(self.read_cigar_file(f)) for f in self.cigar_files]

    def merge_alignments(self, alignments):
        if not alignments:
            return []
            
        merged = []
        current = None
        
        for aln in alignments:
            if aln['op'] == '=' or (not self.show_mismatches and aln['op'] == 'X'):
                if current is None:
                    current = aln.copy()
                elif (current['pos1'] + current['length'] == aln['pos1'] and 
                      current['pos2'] + current['length'] == aln['pos2']):
                    current['length'] += aln['length']
                else:
                    merged.append(current)
                    current = aln.copy()
            else:
                if current is not None:
                    merged.append(current)
                    current = None
                if self.show_mismatches:
                    merged.append(aln)
        
        if current is not None:
            merged.append(current)
            
        return merged

    def plot_region(self, output_file, region_ranges, region_names=None):
        if len(region_ranges) != len(self.regions):
            raise ValueError(f"Number of region ranges ({len(region_ranges)}) must match number of sequences ({len(self.regions)})")
        
        fig_width = 15
        fig, ax = plt.subplots(figsize=(fig_width, 2 * len(self.regions)))
        
        y_positions = np.linspace((len(self.regions) - 1) * 2, 0, len(self.regions))
        rect_height = 0.6
        
        min_x = min(start for start, _ in region_ranges)
        max_x = max(end for _, end in region_ranges)
        
        region_boundaries = []
        for i, region_set in enumerate(self.regions):
            y_pos = y_positions[i]
            region_start, region_end = region_ranges[i]
            
            region_top = y_pos + rect_height/2
            region_bottom = y_pos - rect_height/2
            region_boundaries.append((region_top, region_bottom))
            
            rectangles = [
                Rectangle(
                    (region['start'], y_pos - rect_height/2),
                    region['end'] - region['start'],
                    rect_height,
                    facecolor=boost_rgb(region['rgb'])
                )
                for region in region_set
                if region['end'] > region_start and region['start'] < region_end
            ]
            
            for rect in rectangles:
                ax.add_patch(rect)

        for i, alignment_ops in enumerate(self.alignments):
            _, top_region_bottom = region_boundaries[i]
            lower_region_top, _ = region_boundaries[i + 1]
            
            top_region_start, top_region_end = region_ranges[i]
            bottom_region_start, bottom_region_end = region_ranges[i + 1]

            merged_alignment = self.merge_alignments(alignment_ops)
            
            for aln in merged_alignment:
                if aln['op'] == '=' or (self.show_mismatches and aln['op'] == 'X'):
                    x1_start = aln['pos1']
                    x1_end = aln['pos1'] + aln['length']
                    x2_start = aln['pos2']
                    x2_end = aln['pos2'] + aln['length']
                    
                    if x1_end < top_region_start or x1_start > top_region_end or \
                       x2_end < bottom_region_start or x2_start > bottom_region_end:
                        continue
                    
                    polygon_coords = np.array([
                        [x1_start, top_region_bottom],
                        [x1_end, top_region_bottom],
                        [x2_end, lower_region_top],
                        [x2_start, lower_region_top]
                    ])
                    
                    polygon = plt.Polygon(
                        polygon_coords, 
                        closed=True,
                        facecolor='blue' if aln['op'] == '=' else 'red',
                        alpha=0.2 if aln['op'] == '=' else 0.6,
                        edgecolor='darkred' if aln['op'] == 'X' else None,
                        linewidth=0.5 if aln['op'] == 'X' else 0
                    )
                    ax.add_patch(polygon)

        padding = (max_x - min_x) * 0.05
        ax.set_xlim(min_x - padding, max_x + padding)
        ax.set_ylim(-1, len(self.regions) * 2)
        ax.set_xlabel('Position (bp)')
        ax.set_yticks(y_positions)
        
        if region_names:
            labels = region_names
        else:
            labels = [f"{Path(f).stem} ({start:,}-{end:,})" 
                    for f, (start, end) in zip(self.bed_files, region_ranges)]
        ax.set_yticklabels(labels)
        
        if region_names:
            title = "Zoomed view: " + ", ".join(region_names)
        else:
            title = f"Zoomed region: {min_x:,}-{max_x:,} bp"
        plt.title(title)
        
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        
    def plot_interactive_bokeh(self, html_output):
        from bokeh.plotting import figure, save, output_file as bokeh_output_file
        from bokeh.models import ColumnDataSource, HoverTool, Range1d, CustomJS, Slider
        from bokeh.layouts import column
        
        bokeh_output_file(html_output, title="Interactive Alignment Visualization")
        
        max_length = max(self.sequence_lengths)
        
        track_height = 0.4
        v_spacing = 0.8
        total_height = (len(self.regions) * track_height) + ((len(self.regions) - 1) * v_spacing)
        
        y_positions = []
        track_bottom_positions = []
        connection_points = []
        
        current_y = 0
        
        for i in range(len(self.regions)):
            bottom_y = current_y
            track_bottom_positions.append(bottom_y)
            
            center_y = bottom_y + (track_height / 2)
            y_positions.append(center_y)
            
            if i < len(self.regions) - 1:
                top_connection = bottom_y + track_height
                bottom_connection = bottom_y + track_height + v_spacing
                connection_points.append((top_connection, bottom_connection))
            
            current_y = bottom_y + track_height + v_spacing
        
        p = figure(
            width=1200, 
            height=200 * len(self.regions),
            x_range=(-max_length * 0.02, max_length * 1.02),
            y_range=(0, current_y),
            toolbar_location="above",
            tools="pan,wheel_zoom,box_zoom,reset,save",
            active_scroll="wheel_zoom",
            title="Alignment Visualization"
        )
        
        p.xgrid.visible = False
        p.ygrid.visible = False

        region_hover = HoverTool(
            tooltips=[
                ("Name", "@name"),
                ("Position", "@start - @end"),
                ("Length", "@length bp"),
                ("Color", "@color")
            ],
            renderers=[]
        )
        p.add_tools(region_hover)
        
        # Use shared boost_rgb_hex defined at module level
        
        rect_renderers = []
        
        for i, region_set in enumerate(self.regions):
            bottom_y = track_bottom_positions[i]
            
            rect_data = {
                'x': [region['start'] + (region['end'] - region['start'])/2 for region in region_set],
                'y': [bottom_y + track_height/2 for _ in region_set],
                'width': [region['end'] - region['start'] for region in region_set],
                'height': [track_height for _ in region_set],
                'color': [boost_rgb_hex(region['rgb']) for region in region_set],
                'name': [region['name'] for region in region_set],
                'start': [region['start'] for region in region_set],
                'end': [region['end'] for region in region_set],
                'length': [region['end'] - region['start'] for region in region_set]
            }
            source = ColumnDataSource(data=rect_data)
            
            rect_renderer = p.rect(
                x='x', y='y', width='width', height='height',
                source=source,
                fill_color='color',
                line_color=None,   # no border — cleaner, more vivid colors
                line_width=0
            )
            
            rect_renderers.append(rect_renderer)
            region_hover.renderers.append(rect_renderer)
        
        for i, alignment_ops in enumerate(self.alignments):
            top_y, bottom_y = connection_points[i]
            
            print(f"Track {i} to {i+1} connection: from y={top_y} to y={bottom_y}")
            
            merged_alignment = self.merge_alignments(alignment_ops)
            
            equal_xs = []
            equal_ys = []
            mismatch_xs = []
            mismatch_ys = []
            
            for aln in merged_alignment:
                if aln['op'] == '=' or (self.show_mismatches and aln['op'] == 'X'):
                    x1_start = aln['pos1']
                    x1_end = aln['pos1'] + aln['length']
                    x2_start = aln['pos2']
                    x2_end = aln['pos2'] + aln['length']
                    
                    xs = [x1_start, x1_end, x2_end, x2_start]
                    ys = [top_y, top_y, bottom_y, bottom_y]
                    
                    if aln['op'] == '=':
                        equal_xs.append(xs)
                        equal_ys.append(ys)
                    else:
                        mismatch_xs.append(xs)
                        mismatch_ys.append(ys)
            
            if equal_xs:
                p.patches(
                    xs=equal_xs, 
                    ys=equal_ys, 
                    fill_color='blue', 
                    fill_alpha=0.2, 
                    line_width=0
                )
            
            if mismatch_xs:
                p.patches(
                    xs=mismatch_xs, 
                    ys=mismatch_ys, 
                    fill_color='red', 
                    fill_alpha=0.6, 
                    line_color='darkred',
                    line_width=0.1
                )
        
        # ==========================================
        # ADD DYNAMIC LINE THICKNESS CALLBACK
        # ==========================================
        
        line_callback = CustomJS(args=dict(
            renderers=rect_renderers, 
            x_range=p.x_range,
            max_length=max_length), code="""
            const start = x_range.start;
            const end = x_range.end;
            const view_width = end - start;
            const zoom_factor = view_width / max_length;
            
            let new_line_width;
            let new_line_alpha;
            
            if (zoom_factor > 0.5) {
                new_line_width = 0.1;
                new_line_alpha = 0.2;
            } else if (zoom_factor > 0.2) {
                new_line_width = 0.2;
                new_line_alpha = 0.4;
            } else if (zoom_factor > 0.05) {
                new_line_width = 0.3;
                new_line_alpha = 0.6;
            } else {
                new_line_width = 0.5;
                new_line_alpha = 0.8;
            }
            
            for (let i = 0; i < renderers.length; i++) {
                renderers[i].glyph.line_width = new_line_width;
                renderers[i].glyph.line_alpha = new_line_alpha;
            }
        """)
        
        p.x_range.js_on_change('start', line_callback)
        p.x_range.js_on_change('end', line_callback)
        
        p.xaxis.axis_label = "Position (bp)"
        p.yaxis.ticker = y_positions
        p.yaxis.major_label_overrides = {
            y_pos: f"{Path(f).stem} ({length:,} bp)" 
            for y_pos, f, length in zip(y_positions, self.bed_files, self.sequence_lengths)
        }
        
        x_slider = Slider(
            start=0, 
            end=max_length, 
            value=max_length/2, 
            step=max_length/100,
            title="Horizontal Position (bp)"
        )
        
        slider_callback = CustomJS(args=dict(plot=p, max_length=max_length), code="""
            const value = cb_obj.value;
            const width = plot.x_range.end - plot.x_range.start;
            plot.x_range.start = Math.max(0, value - width/2);
            plot.x_range.end = Math.min(max_length, value + width/2);
        """)
        x_slider.js_on_change('value', slider_callback)
        
        layout = column(p, x_slider)
        save(layout)
        
        print(f"Interactive visualization saved to {html_output}")

    def plot(self, output_file, interactive=False):
        max_length = max(self.sequence_lengths)
        fig_width = min(max(15, max_length / 500), 40)
        fig, ax = plt.subplots(figsize=(fig_width, 2 * len(self.regions)))
        
        y_positions = np.linspace((len(self.regions) - 1) * 2, 0, len(self.regions))
        rect_height = 0.6
        
        region_boundaries = []
        for i, region_set in enumerate(self.regions):
            y_pos = y_positions[i]
            
            region_top = y_pos + rect_height/2
            region_bottom = y_pos - rect_height/2
            region_boundaries.append((region_top, region_bottom))
            
            rectangles = [
                Rectangle(
                    (region['start'], y_pos - rect_height/2),
                    region['end'] - region['start'],
                    rect_height,
                    facecolor=boost_rgb(region['rgb'])
                )
                for region in region_set
            ]

            for rect in rectangles:
                ax.add_patch(rect)

        for i, alignment_ops in enumerate(self.alignments):
            _, top_region_bottom = region_boundaries[i]
            lower_region_top, _ = region_boundaries[i + 1]

            merged_alignment = self.merge_alignments(alignment_ops)
            
            for aln in merged_alignment:
                if aln['op'] == '=' or (self.show_mismatches and aln['op'] == 'X'):
                    x1_start = aln['pos1']
                    x1_end = aln['pos1'] + aln['length']
                    x2_start = aln['pos2']
                    x2_end = aln['pos2'] + aln['length']
                    
                    polygon_coords = np.array([
                        [x1_start, top_region_bottom],
                        [x1_end, top_region_bottom],
                        [x2_end, lower_region_top],
                        [x2_start, lower_region_top]
                    ])
                    
                    polygon = plt.Polygon(
                        polygon_coords, 
                        closed=True,
                        facecolor='blue' if aln['op'] == '=' else 'red',
                        alpha=0.2 if aln['op'] == '=' else 0.6,
                        edgecolor='darkred' if aln['op'] == 'X' else None,
                        linewidth=0.5 if aln['op'] == 'X' else 0
                    )
                    ax.add_patch(polygon)

        ax.set_xlim(-max_length * 0.02, max_length * 1.02)
        ax.set_ylim(-1, len(self.regions) * 2)
        ax.set_xlabel('Position (bp)')
        ax.set_yticks(y_positions)
        
        labels = [f"{Path(f).stem} ({length:,} bp)" 
                  for f, length in zip(self.bed_files, self.sequence_lengths)]
        ax.set_yticklabels(labels)
        
        if interactive:
            plt.subplots_adjust(bottom=0.2)
            
            from matplotlib.widgets import Button, RectangleSelector
            
            def line_select_callback(eclick, erelease):
                x1, y1 = eclick.xdata, eclick.ydata
                x2, y2 = erelease.xdata, erelease.ydata
                ax.set_xlim(min(x1, x2), max(x1, x2))
                ax.set_ylim(min(y1, y2), max(y1, y2))
                fig.canvas.draw_idle()
            
            rs = RectangleSelector(ax, line_select_callback,
                                useblit=True,
                                button=[1],
                                minspanx=5, minspany=5,
                                spancoords='pixels',
                                interactive=True)
            
            ax_reset = plt.axes([0.8, 0.05, 0.1, 0.075])
            button_reset = Button(ax_reset, 'Reset Zoom')
            
            def reset_zoom(event):
                ax.set_xlim(-max_length * 0.02, max_length * 1.02)
                ax.set_ylim(-1, len(self.regions) * 2)
                fig.canvas.draw_idle()
            
            button_reset.on_clicked(reset_zoom)
            
            plt.tight_layout()
            plt.show()
        else:
            plt.savefig(output_file, bbox_inches='tight', dpi=300)
            plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create synteny-style plot from bed files and cigar strings')
    parser.add_argument('--beds', nargs='+', required=True, help='Bed files in order')
    parser.add_argument('--cigars', nargs='+', required=True, help='Cigar string files for pairwise alignments')
    parser.add_argument('--output', required=True, help='Output file name')
    parser.add_argument('--show-mismatches', action='store_true', default=False, help='Show mismatches in the visualization')
    
    output_group = parser.add_argument_group('Output Type')
    output_type = output_group.add_mutually_exclusive_group()
    output_type.add_argument('--static', action='store_true', default=True, help='Create a static image (default)')
    output_type.add_argument('--interactive', action='store_true', help='Show interactive matplotlib plot with zoom capability')
    output_type.add_argument('--web', action='store_true', help='Create interactive web-based visualization with Bokeh')
    
    zoom_group = parser.add_argument_group('Static Zoom Options')
    zoom_group.add_argument('--zoom', action='store_true', help='Create a zoomed-in view of specific regions')
    zoom_group.add_argument('--zoom-regions', nargs='+', help='Comma-separated start,end pairs for each sequence (e.g., "1000,4000")')
    zoom_group.add_argument('--zoom-output', help='Output file for zoomed plot (defaults to adding _zoom suffix)')
    zoom_group.add_argument('--zoom-names', nargs='+', help='Optional names for the zoomed regions')
    
    args = parser.parse_args()
    
    if len(args.cigars) != len(args.beds) - 1:
        raise ValueError("Number of cigar files must be one less than number of bed files")
    
    visualizer = AlignmentVisualizer(args.beds, args.cigars, args.show_mismatches)
    visualizer.process_files()

    if args.web:
        if not args.output.endswith('.html'):
            web_output = args.output + '.html'
        else:
            web_output = args.output
        visualizer.plot_interactive_bokeh(web_output)
    else:
        visualizer.plot(args.output, interactive=args.interactive)
        
        if args.static and args.zoom and args.zoom_regions:
            if len(args.zoom_regions) != len(args.beds):
                raise ValueError(f"Number of zoom regions ({len(args.zoom_regions)}) must match number of bed files ({len(args.beds)})")
            
            region_ranges = []
            for region_str in args.zoom_regions:
                start, end = map(int, region_str.split(','))
                region_ranges.append((start, end))
            
            zoom_output = args.zoom_output if args.zoom_output else args.output.replace('.', '_zoom.')
            visualizer.plot_region(zoom_output, region_ranges, args.zoom_names)

if __name__ == "__main__":
    main()
