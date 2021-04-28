use handlegraph::handle::Handle;
use serde::{Serialize, Deserialize, Serializer, Deserializer};
use std::ops::{Deref, DerefMut};
use serde::de::{self, Visitor};
use std::fmt;
use std::collections::HashMap;
use boomphf::hashmap::NoKeyBoomHashMap;
use boomphf::Mphf;
use crate::kmer::KmerPos;
use std::marker::PhantomData;
use bitvector::*;
use std::fs::File;
use std::io::Write;

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
pub(crate) struct SerializableHandle(u64);

pub fn serialize_object_to_file<T>(object_to_serialize :&T, fileName :&String) -> std::io::Result<()> {
    let serialized_object = bincode::serialize(&object_to_serialize).unwrap();
    let mut file = File::create(fileName)?;
    file.write_all(&serialized_object)?;
    Ok(())
}