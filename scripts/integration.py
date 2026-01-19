# integrate_fastq_support.py
import sys

def integrate_fastq_module():
    """Automatically integrate FASTQ support into improved pipeline."""
    
    # Read the original improved pipeline
    with open('scripts/run_entire_pipeline.py', 'r') as f:
        original_code = f.read()
    
    # Read the FASTQ module
    with open('scripts/fastq_integration_module.py', 'r') as f:
        fastq_module = f.read()
    
    # Find insertion point (after imports, before PART 1)
    insertion_marker = "# ============================================================================\n# PART 1:"
    
    if insertion_marker not in original_code:
        print("ERROR: Could not find insertion point!")
        return False
    
    # Split and insert
    parts = original_code.split(insertion_marker)
    new_code = parts[0] + fastq_module + "\n\n" + insertion_marker + parts[1]
    
    # Write hybrid version
    with open('scripts/run_entire_pipeline_hybrid.py', 'w') as f:
        f.write(new_code)
    
    print("✓ FASTQ module integrated successfully!")
    print("✓ New file created: run_entire_pipeline_hybrid.py")
    return True

if __name__ == "__main__":
    integrate_fastq_module()