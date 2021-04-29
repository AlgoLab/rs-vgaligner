use std::fs::File;
use std::io::Write;

use handlegraph::handle::Handle;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
// Use this struct as a template for serializing the Handle, as Handle does not support Serde
// natively (i.e. the annotation is missing)
pub(crate) struct SerializableHandle(u64);

/// Serialize an object of any kind (as long as it implements Serialize as required by Serde)
/// and put the result of the serialization in a file
pub fn serialize_object_to_file<T: Serialize>(
    object_to_serialize: &T,
    fileName: String,
) -> std::io::Result<()> {
    let serialized_object = bincode::serialize(&object_to_serialize).unwrap();
    let mut file = File::create(fileName)?;
    file.write_all(&serialized_object)?;
    Ok(())
}
