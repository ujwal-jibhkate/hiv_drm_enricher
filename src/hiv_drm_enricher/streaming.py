import pysam
import time
from typing import Iterator

def stream_bam_as_reads(bam_path: str, delay_ms: int = 10) -> Iterator[pysam.AlignedSegment]:
    """
    Opens a BAM file and yields reads one by one to simulate a real-time stream.

    This function acts as a generator, reading a single record at a time,
    making it memory-efficient for very large BAM files. An artificial delay
    is introduced to mimic the pace of a live sequencing run.

    Args:
        bam_path (str): The full path to the input BAM file.
        delay_ms (int): The artificial delay in milliseconds between yielding reads.
                        Set to 0 for no delay.

    Yields:
        Iterator[pysam.AlignedSegment]: An iterator of Pysam AlignedSegment objects.

    Raises:
        FileNotFoundError: If the BAM file or its index cannot be found.
    """
    print(f"INFO: Attempting to open BAM file: {bam_path}")
    try:
        # 'rb' specifies reading a BAM file. Pysam automatically looks for the .bai index.
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            print("INFO: BAM file opened successfully. Starting stream...")
            for read in bamfile:
                # Yield the read object to the caller
                yield read
                
                # Introduce a small delay to simulate real-time arrival
                if delay_ms > 0:
                    time.sleep(delay_ms / 1000.0)

    except FileNotFoundError:
        print(f"ERROR: BAM file not found at '{bam_path}'. Please check the path.")
        raise
    except ValueError as e:
        # Pysam raises ValueError if the index is missing
        print(f"ERROR: Could not open BAM file. Is it indexed? Pysam error: {e}")
        raise