import argparse
import matplotlib.pyplot as plt
import pandas as pd

def parse_region(region):
    """
    Parse region string chr:start-end
    """
    chrom, rest = region.split(":")
    start, end = map(int, rest.split("-"))
    return chrom, start, end

def load_bed(bed_file, chrom, start, end):
    """
    Load BED file and filter intervals overlapping region
    """
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(
        bed_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        names=cols
    )

    df = df[df["chrom"] == chrom]
    df = df[(df["end"] > start) & (df["start"] < end)]
    return df

def plot_tracks(region, bed_files, labels, output):
    chrom, region_start, region_end = parse_region(region)

    fig, ax = plt.subplots(
        figsize=(10, 1.2 * len(bed_files))
    )

    y_height = 0.8

    for i, bed_file in enumerate(bed_files):
        df = load_bed(bed_file, chrom, region_start, region_end)
        y = len(bed_files) - i  # stack top to bottom

        for _, row in df.iterrows():
            ax.add_patch(
                plt.Rectangle(
                    (row["start"], y - y_height / 2),
                    row["end"] - row["start"],
                    y_height,
                    color="black"
                )
            )

        ax.text(
            region_start - (region_end - region_start) * 0.01,
            y,
            labels[i],
            ha="right",
            va="center"
        )

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, len(bed_files) + 1)
    ax.set_yticks([])
    ax.set_xlabel(f"{chrom} position")
    ax.set_title(region)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Plot stacked BED tracks for a genomic region"
    )
    parser.add_argument(
        "--region",
        required=True,
        help="Genomic region (chr:start-end)"
    )
    parser.add_argument(
        "--beds",
        required=True,
        nargs="+",
        help="BED files"
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        help="Track labels (optional)"
    )
    parser.add_argument(
        "--output",
        default="bed_tracks.png",
        help="Output image"
    )

    args = parser.parse_args()

    labels = args.labels if args.labels else args.beds

    plot_tracks(args.region, args.beds, labels, args.output)

if __name__ == "__main__":
    main()
