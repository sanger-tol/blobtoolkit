import groovy.json.JsonSlurper
import java.util.zip.GZIPInputStream

class JSONUtils {
    static String sanitizeBTKMetaJSON(File file) {
        def meta_json = new JsonSlurper().parse(
            new GZIPInputStream(new FileInputStream(file))
        )
        // Remove the Nextflow version, since it will vary between test environments
        meta_json["settings"]["software_versions"].remove("Nextflow")
        // Remove paths to working directories (files created by CAT_CAT)
        meta_json["reads"]["paired"].each { entry ->
            if (entry["file"]?.contains('/work/')) {
                entry.remove('file')
            }
        }
        // Serialise back to a JSON and checksum
        return file.name + ":md5," + groovy.json.JsonOutput.toJson(new TreeMap(meta_json)).md5().toString()
    }
}
