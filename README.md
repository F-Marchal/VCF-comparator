WORK IN PROGRESS README (Files works)

# VCF-comparator
Provide a tool to compare positions presents inside variant call format (VCF files).
This tool provide two way of analyzing VCF files :
- By comparing files in there entirety (Variants summarization).
- By summarizing which variants of each file are the most shared (Files comparison).

## Variants summarization:
### Purpose and functioning
The purpose of this mod is to sumarize the variation contained by a group of files. This is the default mod.

When this mod is used the program will seek each variant of each position of each file. Each position related to a variant are displayed along with a score that represent the proportion of file that contain a variant at this position.

## Files comparison :
### Purpose and functioning
The purpose of this mod is to make the identification of aberrant file easier.  

You can activate this mod with the option `-e` when you call `main.sh`.

When this mod is activated, the program will seek each combinations of files and look how many variants inside the first file match with varaints inside the second file. Each file have two main score attributed to them:
- A percentage of similarity (GSCORE). This represent to what extent the variants inside this files end up in other files and to what extent variants of other files end up in this file.
- A score of inclusion (ISCORE). This score is also percentage. This percentage represent to what extent variants found inside this file can be found inside others files.

At the end of the program a the mean of all GSCORE files is calculated. Furthermore, since each combinations of files is compared a local ISCORE and a local GSORE is calculated for each couple of files.

### Exploiting result
- When a file have ...

