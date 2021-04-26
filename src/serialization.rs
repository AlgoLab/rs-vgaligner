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

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
pub(crate) struct SerializableHandle(u64);

/*
#[derive(Serialize, Deserialize)]
#[serde(remote = "Mphf")]
pub(crate) struct SerializableMphf<T> {
    pub bitvecs: Vec<BitVector>,
    pub ranks: Vec<Vec<u64>>,
    pub phantom: PhantomData<T>,
}

#[derive(Serialize, Deserialize)]
#[serde(remote = "NoKeyBoomHashMap")]
pub(crate) struct SerializableNoKeyBoomHashMap<K,D1> {
    pub mphf: SerializableMphf<K>,
    pub values: Vec<D1>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct KmerTable {
    #[serde(with = "SerializableNoKeyBoomHashMap")]
    pub table: NoKeyBoomHashMap<u64, KmerPos>,
    pub size: u64
}
 */