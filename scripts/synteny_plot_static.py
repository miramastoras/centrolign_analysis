#!/usr/bin/env python3

import argparse
import colorsys
import re
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.transforms import blended_transform_factory
import numpy as np


def boost_rgb(rgb_str, saturation_boost=2.0, value_floor=0.85):
    try:
        r, g, b = [int(x) / 255.0 for x in rgb_str.split(',')]
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        s = min(1.0, s * saturation_boost)
        v = max(value_floor, v)
        return colorsys.hsv_to_rgb(h, s, v)
    except Exception:
        return (0.8, 0.8, 0.8)


def boost_rgb_colorpop(rgb_str):
    """Boost colors vividly, but preserve earthy/muted tones (browns, tans) as distinguishable."""
    try:
        r, g, b = [int(x) / 255.0 for x in rgb_str.split(',')]
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        if s < 0.55 and v < 0.70:
            # Earthy/muted tone — boost gently to preserve character
            s = min(0.75, s * 1.5)
            v = min(0.75, v * 1.3)
        else:
            # Vivid or bright color — push to neon
            s = min(1.0, s * 3.0)
            v = max(0.92, v)
        return colorsys.hsv_to_rgb(h, s, v)
    except Exception:
        return (0.8, 0.8, 0.8)


class AlignmentVisualizer:
    def __init__(self, bed_files, cigar_files, show_mismatches, colorpop=False):
        self.bed_files = bed_files
        self.cigar_files = cigar_files
        self.regions = []
        self.alignments = []
        self.show_mismatches = show_mismatches
        self.colorpop = colorpop
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
                    'name': fields[3] if len(fields) > 3 else '.',
                    'score': float(fields[4]) if len(fields) > 4 else 0.0,
                    'strand': fields[5] if len(fields) > 5 else '.',
                    'thick_start': int(fields[6]) if len(fields) > 6 else int(fields[1]),
                    'thick_end': int(fields[7]) if len(fields) > 7 else int(fields[2]),
                    'rgb': fields[8] if len(fields) > 8 else '128,128,128'
                }
                for line in lines
                if not line.startswith('#') and len((fields := line.strip().split('\t'))) >= 3
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
            self.sequence_lengths.append(
                max(region['end'] for region in regions) if regions else 0
            )

        self.alignments = [
            self.parse_cigar_string(self.read_cigar_file(f))
            for f in self.cigar_files
        ]

    def _cigar_track_limits(self):
        """
        Return the sequence length consumed by each track as defined by adjacent cigars.
        For each cigar, compute final pos1 (seq1 length) and pos2 (seq2 length).
        Track 0       -> cigar[0] pos1
        Track i (mid) -> max(cigar[i-1] pos2, cigar[i] pos1)
        Track n-1     -> cigar[-1] pos2
        This naturally handles cigars ending in gaps: the gapped sequence still
        gets its full cigar-implied length.
        """
        cigar_lens = []
        for ops in self.alignments:
            if not ops:
                cigar_lens.append((0, 0))
                continue
            last = ops[-1]
            pos1 = last['pos1'] + (last['length'] if last['op'] in ('=', 'X', 'D') else 0)
            pos2 = last['pos2'] + (last['length'] if last['op'] in ('=', 'X', 'I') else 0)
            cigar_lens.append((pos1, pos2))

        n = len(self.regions)
        limits = []
        for i in range(n):
            if i == 0:
                limits.append(cigar_lens[0][0])
            elif i == n - 1:
                limits.append(cigar_lens[-1][1])
            else:
                limits.append(max(cigar_lens[i - 1][1], cigar_lens[i][0]))
        return limits

    @staticmethod
    def _short_label(name):
        """Extract label before first '(' — e.g. 'gSat(GSATII,TAR1)' -> 'gSat'."""
        return name.split('(')[0].strip()

    def _legend_entries(self):
        """Return ordered list of (rgb_tuple, label) with no duplicates by rgb."""
        seen = {}  # rgb_str -> label (first occurrence wins)
        for region_set in self.regions:
            for region in region_set:
                rgb_str = region['rgb']
                if rgb_str not in seen:
                    seen[rgb_str] = self._short_label(region['name'])
        # Build handles sorted by label name
        entries = sorted(seen.items(), key=lambda x: x[1].lower())
        if self.colorpop:
            return [(boost_rgb_colorpop(rgb), label) for rgb, label in entries]
        return [
            (tuple(int(v) / 255 for v in rgb.split(',')), label)
            for rgb, label in entries
        ]

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

    def plot(self, output_file, track_labels=None, alignment_labels=None):
        n_tracks = len(self.regions)
        track_limits = self._cigar_track_limits()
        max_length = max(track_limits)

        # Square-ish figure: width based on sequence length, height per track
        track_height_in = 1.5   # inches per track
        gap_height_in = 1.2     # inches per gap between tracks
        total_height = n_tracks * track_height_in + (n_tracks - 1) * gap_height_in + 1.0

        # Scale fonts and legend elements relative to a 3-track reference
        ref_height = 3 * track_height_in + 2 * gap_height_in + 1.0  # 7.9 in
        scale = total_height / ref_height
        fontsize = max(8, round(22 * scale))
        handle_size = round(1.5 * scale * 10) / 10
        handle_height = round(0.8 * scale * 10) / 10  # keeps squares from overlapping rows

        fig_width = min(max(total_height * 2.7, 16), 60)  # wider than tall
        fig, ax = plt.subplots(figsize=(fig_width, total_height))

        # Y layout: tracks separated by gaps, bottom track at y=0
        # Each track occupies [y_center - rect_h/2, y_center + rect_h/2]
        unit = track_height_in + gap_height_in  # y-units per track
        rect_height = 0.6 * track_height_in / unit  # scale rect height to data coords

        # y positions in data coords (top track highest y value)
        y_positions = [(n_tracks - 1 - i) * 2.0 for i in range(n_tracks)]
        rect_height = 0.7  # in data units

        # Store exact region boundaries
        region_boundaries = []
        for i, region_set in enumerate(self.regions):
            y_pos = y_positions[i]
            region_top = y_pos + rect_height / 2
            region_bottom = y_pos - rect_height / 2
            region_boundaries.append((region_top, region_bottom))
            limit = track_limits[i]

            for region in region_set:
                if region['start'] >= limit:
                    continue
                clipped_end = min(region['end'], limit)
                color = boost_rgb_colorpop(region['rgb']) if self.colorpop else [int(x) / 255 for x in region['rgb'].split(',')]
                rect = plt.Rectangle(
                    (region['start'], region_bottom),
                    clipped_end - region['start'],
                    rect_height,
                    facecolor=color,
                    edgecolor='none',
                    linewidth=0
                )
                ax.add_patch(rect)

        # Plot alignments as polygons between adjacent tracks
        for i, alignment_ops in enumerate(self.alignments):
            _, top_region_bottom = region_boundaries[i]
            lower_region_top, _ = region_boundaries[i + 1]

            if alignment_labels and i < len(alignment_labels):
                gap_mid_y = (top_region_bottom + lower_region_top) / 2
                trans = blended_transform_factory(ax.transAxes, ax.transData)
                ax.text(
                    -0.01, gap_mid_y, alignment_labels[i],
                    transform=trans,
                    ha='right', va='center',
                    fontsize=round(fontsize * 1.25), fontweight='bold',
                    color='#444444', zorder=5
                )

            merged = self.merge_alignments(alignment_ops)

            for aln in merged:
                if aln['op'] == '=' or (self.show_mismatches and aln['op'] == 'X'):
                    x1_start = aln['pos1']
                    x1_end = aln['pos1'] + aln['length']
                    x2_start = aln['pos2']
                    x2_end = aln['pos2'] + aln['length']

                    polygon_coords = np.array([
                        [x1_start, top_region_bottom],
                        [x1_end,   top_region_bottom],
                        [x2_end,   lower_region_top],
                        [x2_start, lower_region_top]
                    ])

                    polygon = plt.Polygon(
                        polygon_coords,
                        closed=True,
                        facecolor='blue' if aln['op'] == '=' else 'red',
                        alpha=0.2 if aln['op'] == '=' else 0.6,
                        edgecolor='darkred' if aln['op'] == 'X' else None,
                        linewidth=0.5 if aln['op'] == 'X' else 0,
                        rasterized=True
                    )
                    ax.add_patch(polygon)

        # Axes limits
        padding = max_length * 0.01
        ax.set_xlim(-padding, max_length + padding)
        y_min = y_positions[-1] - rect_height / 2 - 0.3
        y_max = y_positions[0] + rect_height / 2 + 0.3
        ax.set_ylim(y_min, y_max)

        # Y-axis track labels
        ax.set_yticks(y_positions)
        if track_labels and len(track_labels) == n_tracks:
            labels = track_labels
        else:
            labels = [
                f"{Path(f).stem}\n({length:,} bp)"
                for f, length in zip(self.bed_files, self.sequence_lengths)
            ]
        ax.set_yticklabels(labels, fontsize=fontsize, fontweight='bold')

        # X-axis formatting
        ax.set_xlabel('Position (bp)', fontsize=fontsize, labelpad=10)
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f} Mb' if max_length > 1e6 else f'{x:,.0f}')
        )
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize, length=0)

        # Remove top/right spines for cleaner look
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Annotation legend (top-right)
        legend_entries = self._legend_entries()
        annot_handles = [
            mpatches.Patch(facecolor=rgb, edgecolor='none', label=label)
            for rgb, label in legend_entries
        ]
        leg1 = ax.legend(
            handles=annot_handles,
            loc='upper left',
            bbox_to_anchor=(1.01, 1.0),
            fontsize=fontsize,
            frameon=False,
            handlelength=handle_size,
            handleheight=handle_height,
            borderpad=1.0,
            labelspacing=1.0,
            title='Annotation',
            title_fontsize=fontsize
        )
        ax.add_artist(leg1)
        leg1.get_title().set_fontweight('bold')

        # Alignment legend (below annotation legend)
        match_color = (0.8, 0.8, 1.0)  # blue @ alpha=0.2 on white
        match_label = 'match' if self.show_mismatches else 'aligned'
        aln_handles = [
            mpatches.Patch(facecolor='white', edgecolor='black', linewidth=1.5, label='gap'),
            mpatches.Patch(facecolor=match_color, edgecolor='none', label=match_label),
        ]
        if self.show_mismatches:
            aln_handles.append(
                mpatches.Patch(facecolor='darkred', edgecolor='darkred', linewidth=0.5, label='mismatch')
            )
        ax.legend(
            handles=aln_handles,
            loc='lower left',
            bbox_to_anchor=(1.01, 0.0),
            fontsize=fontsize,
            frameon=False,
            handlelength=handle_size,
            handleheight=round(handle_height * 0.5 * 10) / 10,
            borderpad=1.0,
            labelspacing=1.0,
            title='Alignment',
            title_fontsize=fontsize
        )
        ax.get_legend().get_title().set_fontweight('bold')


        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved to {output_file}")

    def plot_region(self, output_file, region_ranges, track_labels=None):
        """Plot a zoomed region for each track."""
        if len(region_ranges) != len(self.regions):
            raise ValueError(
                f"Number of region ranges ({len(region_ranges)}) must match "
                f"number of sequences ({len(self.regions)})"
            )

        n_tracks = len(self.regions)
        fig_width = 12
        total_height = n_tracks * 1.5 + (n_tracks - 1) * 1.2 + 1.0
        fig, ax = plt.subplots(figsize=(fig_width, total_height))

        y_positions = [(n_tracks - 1 - i) * 1.0 for i in range(n_tracks)]
        rect_height = 0.55

        min_x = min(s for s, _ in region_ranges)
        max_x = max(e for _, e in region_ranges)

        region_boundaries = []
        for i, region_set in enumerate(self.regions):
            y_pos = y_positions[i]
            region_start, region_end = region_ranges[i]
            region_top = y_pos + rect_height / 2
            region_bottom = y_pos - rect_height / 2
            region_boundaries.append((region_top, region_bottom))

            for region in region_set:
                if region['end'] <= region_start or region['start'] >= region_end:
                    continue
                color = boost_rgb_colorpop(region['rgb']) if self.colorpop else [int(x) / 255 for x in region['rgb'].split(',')]
                rect = plt.Rectangle(
                    (region['start'], region_bottom),
                    region['end'] - region['start'],
                    rect_height,
                    facecolor=color,
                    edgecolor='none',
                    linewidth=0
                )
                ax.add_patch(rect)

        for i, alignment_ops in enumerate(self.alignments):
            _, top_region_bottom = region_boundaries[i]
            lower_region_top, _ = region_boundaries[i + 1]
            top_start, top_end = region_ranges[i]
            bot_start, bot_end = region_ranges[i + 1]

            merged = self.merge_alignments(alignment_ops)
            for aln in merged:
                if aln['op'] == '=' or (self.show_mismatches and aln['op'] == 'X'):
                    x1_s = aln['pos1']
                    x1_e = aln['pos1'] + aln['length']
                    x2_s = aln['pos2']
                    x2_e = aln['pos2'] + aln['length']

                    if x1_e < top_start or x1_s > top_end or \
                       x2_e < bot_start or x2_s > bot_end:
                        continue

                    polygon_coords = np.array([
                        [x1_s, top_region_bottom],
                        [x1_e, top_region_bottom],
                        [x2_e, lower_region_top],
                        [x2_s, lower_region_top]
                    ])
                    polygon = plt.Polygon(
                        polygon_coords, closed=True,
                        facecolor='blue' if aln['op'] == '=' else 'red',
                        alpha=0.2 if aln['op'] == '=' else 0.6,
                        edgecolor='darkred' if aln['op'] == 'X' else None,
                        linewidth=0.5 if aln['op'] == 'X' else 0
                    )
                    ax.add_patch(polygon)

        padding = (max_x - min_x) * 0.02
        ax.set_xlim(min_x - padding, max_x + padding)
        ax.set_ylim(
            y_positions[-1] - rect_height / 2 - 0.3,
            y_positions[0] + rect_height / 2 + 0.3
        )
        ax.set_xlabel('Position (bp)', fontsize=14, labelpad=8)
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f'{x/1e6:.2f} Mb' if (max_x - min_x) > 1e5 else f'{x:,.0f}')
        )
        ax.tick_params(axis='x', labelsize=12)
        ax.set_yticks(y_positions)
        if track_labels and len(track_labels) == n_tracks:
            labels = track_labels
        else:
            labels = [
                f"{Path(f).stem}\n({s:,}–{e:,})"
                for f, (s, e) in zip(self.bed_files, region_ranges)
            ]
        ax.set_yticklabels(labels, fontsize=14, fontweight='bold')
        ax.tick_params(axis='y', length=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved zoomed plot to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create synteny-style plot from BED files and CIGAR strings (static matplotlib)'
    )
    parser.add_argument('--beds', nargs='+', required=True,
                        help='BED files in order (top to bottom)')
    parser.add_argument('--cigars', nargs='+', required=True,
                        help='CIGAR string files for pairwise alignments between adjacent tracks')
    parser.add_argument('--output', required=True,
                        help='Output file name (e.g. plot.pdf or plot.png)')
    parser.add_argument('--show-mismatches', action='store_true', default=False,
                        help='Highlight mismatches (X operations) in red')
    parser.add_argument('--colorpop', action='store_true', default=False,
                        help='Boost BED track colors to be more vivid and neon')
    parser.add_argument('--labels', nargs='+', default=None,
                        help='Custom labels for each BED track (top to bottom). '
                             'If not given, uses file stem + sequence length.')
    parser.add_argument('--alignment-labels', nargs='+', default=None,
                        help='Labels for each alignment gap between tracks (top to bottom). '
                             'Must have len(beds)-1 entries.')

    # Zoom options
    zoom_group = parser.add_argument_group('Zoom Options')
    zoom_group.add_argument('--zoom', action='store_true',
                            help='Also create a zoomed-in view of specific regions')
    zoom_group.add_argument('--zoom-regions', nargs='+',
                            help='Comma-separated start,end pairs for each track '
                                 '(e.g. "1000,4000" "2000,5000")')
    zoom_group.add_argument('--zoom-output',
                            help='Output file for zoomed plot (default: adds _zoom suffix)')

    args = parser.parse_args()

    if len(args.cigars) != len(args.beds) - 1:
        parser.error(
            f"Need exactly len(beds)-1 cigar files. "
            f"Got {len(args.beds)} beds and {len(args.cigars)} cigars."
        )

    if args.labels and len(args.labels) != len(args.beds):
        parser.error(
            f"--labels must have the same number of entries as --beds "
            f"({len(args.beds)} expected, got {len(args.labels)})."
        )

    if args.alignment_labels and len(args.alignment_labels) != len(args.beds) - 1:
        parser.error(
            f"--alignment-labels must have len(beds)-1 entries "
            f"({len(args.beds) - 1} expected, got {len(args.alignment_labels)})."
        )

    visualizer = AlignmentVisualizer(args.beds, args.cigars, args.show_mismatches, colorpop=args.colorpop)
    visualizer.process_files()

    visualizer.plot(args.output, track_labels=args.labels, alignment_labels=args.alignment_labels)

    if args.zoom and args.zoom_regions:
        if len(args.zoom_regions) != len(args.beds):
            parser.error(
                f"--zoom-regions must have one entry per bed file "
                f"({len(args.beds)} expected, got {len(args.zoom_regions)})."
            )
        region_ranges = [tuple(map(int, r.split(','))) for r in args.zoom_regions]
        zoom_out = args.zoom_output or args.output.rsplit('.', 1)[0] + '_zoom.' + args.output.rsplit('.', 1)[-1]
        visualizer.plot_region(zoom_out, region_ranges, track_labels=args.labels)


if __name__ == "__main__":
    main()
