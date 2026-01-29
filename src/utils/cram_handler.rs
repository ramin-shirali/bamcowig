use noodles_cram as cram;

fn cramReader() {
    let mut reader = cram::io::reader::Builder::default().build_from_path("sample.bam")?;
    let header = reader.read_header()?;

    for result in reader.records() {
        let record = result?; 
        // ...
    }
    // Ok::<_, io::Error>(())
}