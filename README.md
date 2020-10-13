# BSBoltManuscript
Supporting Material and Analysis Code for *BiSulfite Bolt: A BiSulfite Sequencing Analysis Platform*

## Manuscript Simulation Data

The simulation data used in the mansuscript is deposited in a publically accessible Google Drive, [BSBoltSimulationData](https://drive.google.com/drive/u/1/folders/16bHclPwS2_-7cheRuT8hJKOBZAsHse-f). The simulation conditions used and the output directory are provided in the table below. Each simulation directory contains compressed fastq file(s) and serialized python objects reporting the methylation status of all methylatable bases within each simulation contig.

|   Average Read Depth |   Mutation Rate |   Sequencing Error | Sequencing Type   | Library Type   | Simulation Directory            |
|---------------------:|----------------:|-------------------:|:------------------|:---------------|:--------------------------------|
|                   20 |           0.005 |              0.005 | Paired End        | Directional    | pe_directional_50               |
|                   20 |           0.005 |              0.005 | Paired End        | Directional    | pe_directional_100              |
|                   20 |           0.005 |              0.005 | Paired End        | Directional    | pe_directional_150              |
|                   30 |           0.005 |              0.005 | Paired End        | Undirectional  | pe_undirectional_50             |
|                   30 |           0.005 |              0.005 | Paired End        | Undirectional  | pe_undirectional_100            |
|                   30 |           0.005 |              0.005 | Paired End        | Undirectional  | pe_undirectional_150            |
|                   20 |           0.005 |              0.005 | Single End        | Directional    | se_directional_50               |
|                   20 |           0.005 |              0.005 | Single End        | Directional    | se_directional_100              |
|                   20 |           0.005 |              0.005 | Single End        | Directional    | se_directional_150              |
|                   30 |           0.005 |              0.005 | Single End        | Undirectional  | se_undirectional_50             |
|                   30 |           0.005 |              0.005 | Single End        | Undirectional  | se_undirectional_100            |
|                   30 |           0.005 |              0.005 | Single End        | Undirectional  | se_undirectional_150            |
|                    8 |           0.005 |              0.005 | Paired End        | Directional    | pe_low_coverage_directional_50  |
|                    8 |           0.005 |              0.005 | Paired End        | Directional    | pe_low_coverage_directional_100 |
|                    8 |           0.005 |              0.005 | Paired End        | Directional    | pe_low_coverage_directional_150 |
|                    8 |           0.005 |              0.005 | Single End        | Directional    | se_low_coverage_directional_50  |
|                    8 |           0.005 |              0.005 | Single End        | Directional    | se_low_coverage_directional_100 |
|                    8 |           0.005 |              0.005 | Single End        | Directional    | se_low_coverage_directional_150 |
|                    8 |           0.01  |              0.02  | Paired End        | Directional    | pe_error_directional_50         |
|                    8 |           0.01  |              0.02  | Paired End        | Directional    | pe_error_directional_100        |
|                    8 |           0.01  |              0.02  | Paired End        | Directional    | pe_error_directional_150        |