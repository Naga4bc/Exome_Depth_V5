# ExomeDepth CNV Analysis Pipeline

A dockerized R pipeline for detecting Copy Number Variations (CNVs) from exome sequencing data using the ExomeDepth algorithm.

## Overview

This pipeline automates CNV detection from BAM files using ExomeDepth, with support for multiple sample comparisons and gene-specific analysis. It's designed to run in a Docker container for reproducibility and ease of deployment.

## Features

- **Automated Docker deployment**: Self-launching Docker container with proper user permissions
- **Flexible sample processing**: Supports multiple samples with comparison controls
- **Multi-platform BED file support**: Auto-selects appropriate target regions based on sample naming
- **Gene-specific analysis**: Generates plots and analysis for specified genes
- **Comprehensive logging**: Timestamped logs for monitoring and debugging
- **Batch processing**: Processes multiple sample sets from a single input file

## Requirements

### System Requirements
- Linux system with Docker installed
- User with sudo privileges for Docker execution
- Minimum 8GB RAM recommended
- Storage space for BAM files and output results

### Docker Image
- `euformatics/exomedepth:v1.1` (automatically pulled)

## Installation

1. Clone or download the script:
```bash
# Download the main script
wget [script-url] -O Exome_Depth_CNACall_V5.R
```

2. Ensure Docker is installed and running:
```bash
sudo systemctl start docker
sudo systemctl enable docker
```

3. Set up the required directory structure:
```bash
mkdir -p /home/bioinfo4/Patient_Samples/Exome_Depth_Dockerization
cd /home/bioinfo4/Patient_Samples/Exome_Depth_Dockerization
```

## Input Files

### Required Files

1. **`sample_genelist.txt`** - Sample and gene specifications
   ```
   SAMPLE1,SAMPLE2,SAMPLE3:GENE1,GENE2,GENE3
   SAMPLE4,SAMPLE5,SAMPLE6:GENE4,GENE5
   ```
   Format: `main_sample,control_sample1,control_sample2:gene1,gene2,gene3`

2. **`projects.txt`** - List of BaseSpace project names
   ```
   Project_Name_1
   Project_Name_2
   Project_Name_3
   ```

3. **BAM files** - Located in BaseSpace project structure:
   ```
   /workspace/basespace/Projects/[PROJECT_NAME]/AppResults/[SAMPLE_ID]/Files/[SAMPLE_ID].bam
   ```

4. **BED files** - Target region definitions (auto-selected based on sample naming):
   - `Indiegene_Target_2109PD006-V1_4BaseCare_1K_DNA_GRCh37.bed` (for -CE- or -CEFu- samples)
   - `TarGT_First_v2_CDS_GRCh37_13_Mar_23.bed` (for -FEV2F2both- samples)
   - `T1_CDS_V3_hg19_9D24.bed` (for -CDSV36XFUS- samples)
   - `SureSelect_V8_Coverde_Modified.bed` (for -SE8- samples)

## Usage

### Basic Usage
```bash
# Navigate to the working directory
cd /home/bioinfo4/Patient_Samples/Exome_Depth_Dockerization

# Run the pipeline
Rscript Exome_Depth_CNACall_V5.R
```

The script will automatically:
1. Launch itself in a Docker container if not already running in one
2. Read input files (`sample_genelist.txt` and `projects.txt`)
3. Process each sample set sequentially
4. Generate outputs in organized folder structure

### Sample Naming Conventions

The pipeline auto-selects BED files based on sample ID patterns:
- **-CE-** or **-CEFu-**: Uses Indiegene target regions
- **-FEV2F2both-**: Uses TarGT First v2 CDS regions
- **-CDSV36XFUS-**: Uses T1 CDS V3 regions
- **-SE8-**: Uses SureSelect V8 regions

## Output Structure

```
/workspace/
├── [SAMPLE_ID]/
│   ├── [SAMPLE_ID]_ED_CNA.csv           # CNV calls
│   └── [GENE_NAME]/
│       └── [SAMPLE_ID]_plot_[GENE_NAME]_[START]_[END]_Amplification.png
├── ExomeDepth_log.txt                    # Pipeline logs
├── sample_genelist.txt                   # Input file
└── projects.txt                          # Input file
```

### Output Files

- **`[SAMPLE_ID]_ED_CNA.csv`**: Complete CNV calls with statistics
- **Gene plots**: PNG files showing CNV regions for each specified gene
- **`ExomeDepth_log.txt`**: Timestamped execution log

## Configuration

### Key Parameters
- **Threshold**: Fixed at 0.0001 (transition probability for CNV calling)
- **BED file selection**: Automatic based on sample naming patterns
- **Reference samples**: Uses comparison samples from the same line in `sample_genelist.txt`

### Modifying Parameters
To change the threshold or other parameters, edit the script:
```r
threshold <- 0.0001  # Adjust as needed
```

## Algorithm Details

The pipeline uses the ExomeDepth algorithm which:
1. Counts reads in target regions from BAM files
2. Selects optimal reference samples for comparison
3. Models read count ratios using a hidden Markov model
4. Identifies regions with significant copy number changes
5. Generates statistical confidence measures for each call

## Troubleshooting

### Common Issues

1. **"BAM file not found"**
   - Check that BAM files exist in the expected BaseSpace structure
   - Verify project names in `projects.txt`
   - Ensure sample IDs match exactly

2. **"No matching BED file found"**
   - Check sample naming conventions
   - Ensure required BED files are present in `/workspace/`

3. **Docker permission errors**
   - Verify user has sudo privileges
   - Check Docker service is running: `sudo systemctl status docker`

4. **Memory issues**
   - Increase available RAM
   - Consider processing fewer samples simultaneously

### Log Analysis
Check `ExomeDepth_log.txt` for detailed execution information:
```bash
tail -f /workspace/ExomeDepth_log.txt
```

## Version History

- **V5** (Sep 17, 2025): 
  - Removed CLI arguments
  - Added automatic Docker activation
  - Improved logging and error handling
  - Auto-selection of BED files and parameters

- **Original** (Dec 21, 2024): Initial version by Vyomesh J & Bhumika P

## Dependencies

### R Packages
- `ExomeDepth`: CNV detection algorithm

### System Dependencies
- Docker runtime
- Linux environment with BaseSpace data structure

## Citation

If you use this pipeline in your research, please cite:
- ExomeDepth: Plagnol V, Curtis J, Epstein M, et al. A robust model for read count data in exome sequencing experiments and implications for copy number variant calling. Bioinformatics. 2012;28(21):2747-2754.

## Support

For issues and questions:
1. Check the log file for detailed error messages
2. Verify input file formats match the specifications
3. Ensure all required files are in the correct locations

## License

Please check with the original authors for licensing information.
