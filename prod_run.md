To run sanger-tol/blobtoolkit v0.2.0 in production, follow these steps:
- Format the samplesheet
  ```
  sample,datatype,datafile
  ERR####,hic,/path/to/file
  ERR####,illumina,/path/to/file
  ERR####,ont,/path/to/file
  ERR####,pacbio,/path/to/file
  ```
- You can provide as many samples as you want, the pipeline will select the top 3 read datasets
  - To increase the limit add `--read-runs #` argument for `blobtoolkit/config` module
- Provide the input parameter matching either the `test` profile or `test_raw` profile if running with alignment.
  - This will trigger the pipeline to create the appropriate config and to update with software versions in the end
- Once the pipeline is complete, copy the data to ToL directory structure similar to other pipelines.
- Once the pipeline is complete, tar the blobdir only and add to a central `/lustre` location similar to HiGlass
  - This is best done as part of the Prefect script rather than the pipeline as this will change with the new GAP portal.
  - Rich will then pull the data from the farm, no need to send anything to the production server.
- This data is also meant to be backed on s3 as part of the GAP portal, so please save all files
  - Especially the params.json and execution trace along with the blobdir and BUSCO results.
