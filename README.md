Takes JSON. constructs PATRIC tn-seq run using TRANSIT

```
usage: p3_tnseq.py [-h] (--jfile JFILE | --jstring JSTRING)
                   (--sfile SFILE | --sstring SSTRING) [-o O]

optional arguments:
  -h, --help         show this help message and exit
  --jfile JFILE      json file for job
                     {"reference_genome_id":[x],"experimental_conditions":[x],
                     "transit_params":{key:value}, "output_path":x,
                     "read_files":x
  --jstring JSTRING  json string from user input
  --sfile SFILE      server setup JSON file
  --sstring SSTRING  server setup JSON string
  -o O               output directory. defaults to current directory.
```
