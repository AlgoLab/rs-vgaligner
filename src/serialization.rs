use std::fs::File;
use std::io::{BufReader, Write};

use handlegraph::handle::Handle;
use json::JsonValue;
//use rayon::prelude::*;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[derive(Serialize, Deserialize)]
//#[serde(remote = "Handle")]
// Use this struct as a template for serializing the Handle, as Handle does not support Serde
// natively (i.e. the annotation is missing)
pub(crate) struct SerializableHandle(u64);

// Custom serializer for Index.edges (of type Vec<Handle> which does not implement Serialize)
pub fn handle_vec_ser<S: Serializer>(vec: &Vec<Handle>, serializer: S) -> Result<S::Ok, S::Error> {
    let vec2: Vec<SerializableHandle> = vec
        .iter()
        .map(|h| SerializableHandle(h.as_integer()))
        .collect();
    vec2.serialize(serializer)
}
// Custom deserializer for Index.edges (of type Vec<Handle> which does not implement Deserialize)
pub fn handle_vec_deser<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<Handle>, D::Error> {
    let vec: Vec<SerializableHandle> = Deserialize::deserialize(deserializer)?;
    Ok(vec
        .iter()
        .map(|sh| Handle::from_integer(sh.0))
        .collect())
}

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

/// Store path-to-sequence mappings is a JSON file
pub fn store_mappings_in_file(object: &String, file_name: String) -> std::io::Result<()> {
    let mut file =
        File::create(&file_name).unwrap_or_else(|_| panic!("Couldn't create file {}", &file_name));
    file.write_all(object.as_str().as_ref())
        .unwrap_or_else(|_| panic!("Couldn't write to file {}", &file_name));
    Ok(())
}
