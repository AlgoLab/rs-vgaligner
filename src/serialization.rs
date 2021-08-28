use std::fs::File;
use std::io::{BufReader, Write};

use handlegraph::handle::Handle;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
// Use this struct as a template for serializing the Handle, as Handle does not support Serde
// natively (i.e. the annotation is missing)
pub(crate) struct SerializableHandle(u64);

/// Serialize an object of any kind (as long as it implements Serialize as required by Serde)
/// and put the result of the serialization in a file
pub fn serialize_object_to_file<T: Serialize>(
    object_to_serialize: &T,
    file_name: String,
) -> std::io::Result<()> {
    let serialized_object = bincode::serialize(&object_to_serialize).unwrap();
    let mut file =
        File::create(&file_name).unwrap_or_else(|_| panic!("Couldn't create file {}", &file_name));
    file.write_all(&serialized_object)
        .unwrap_or_else(|_| panic!("Couldn't write to file {}", &file_name));
    Ok(())
}

/// Deserialize an object of any kind (as long as it implements Deserialize as required by Serde)
/// and returns the result of the deserialization
pub fn deserialize_object_from_file<T>(file_name: String) -> T
where
    T: DeserializeOwned,
{
    let reader = BufReader::new(File::open(file_name).unwrap());
    bincode::deserialize_from(reader).unwrap()
}
