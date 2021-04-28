use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};

use bitvector::*;
use boomphf::hashmap::NoKeyBoomHashMap;
use boomphf::Mphf;
use handlegraph::handle::Handle;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde::de::{self, Visitor};

use crate::kmer::KmerPos;

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
pub(crate) struct SerializableHandle(u64);

pub fn serialize_object_to_file<T: Serialize>(
    object_to_serialize: &T,
    fileName: String,
) -> std::io::Result<()> {
    let serialized_object = bincode::serialize(&object_to_serialize).unwrap();
    let mut file = File::create(fileName)?;
    file.write_all(&serialized_object)?;
    Ok(())
}
