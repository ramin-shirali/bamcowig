use noodles_bam as bam;

fn bamReader() {
    let mut reader = bam::io::reader::Builder::default().build_from_path("sample.bam")?;
    let header = reader.read_header()?;

    for result in reader.records() {
        let record = result?;
        // ...
    }
    // Ok::<_, io::Error>(())
}