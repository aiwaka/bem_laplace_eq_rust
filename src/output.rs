use crate::type_and_trait::Point;
use std::fs::File;
use std::io::Write;

const OUTPUT_FILE_NAME: &str = "plot.txt";

pub fn output_data(data_list: &[(Point<f64>, f64)]) -> Result<(), Box<dyn std::error::Error>> {
    let mut fw = File::create(OUTPUT_FILE_NAME).unwrap();
    for data in data_list.iter() {
        #[allow(clippy::write_with_newline)]
        write!(fw, "{:.04}\t{:.04}\t{:.04}\n", data.0[0], data.0[1], data.1)?;
    }
    Ok(())
}
