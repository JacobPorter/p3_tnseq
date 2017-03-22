Takes JSON. constructs PATRIC tn-seq run using TRANSIT

```
usage: p3_tnseq.py [-h] (--ufile UFILE | --ustring USTRING)
                   (--sfile SFILE | --sstring SSTRING) [-o O]

optional arguments:
  -h, --help         show this help message and exit
  --ufile UFILE      json file for job
                     {"reference_genome_id":[x],"experimental_conditions":[x],
                     "transit_params":{key:value}, "output_path":x,
                     "read_files":x
  --ustring USTRING  json string from user input
  --sfile SFILE      server setup JSON file
  --sstring SSTRING  server setup JSON string
  -o O               output directory. defaults to current directory.
```
