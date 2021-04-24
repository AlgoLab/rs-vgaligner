use handlegraph::handle::Handle;
use serde::{Serialize, Deserialize, Serializer, Deserializer};
use std::ops::{Deref, DerefMut};
use serde::de::{self, Visitor};
use std::fmt;

/// Create a wrapper for Handle so that it can be serialized
#[derive(Debug, Clone)]
pub(crate) struct SerializableHandle(Handle);

#[derive(Serialize, Deserialize)]
#[serde(remote = "Handle")]
pub(crate) struct SerialHandle(u64);

/// Serialization trait
impl Serialize for SerializableHandle {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
    {
        serializer.serialize_u64(self.0.0)
    }
}

/// Deserialization trait
struct SerializableHandleVisitor;

impl<'de> Visitor<'de> for SerializableHandleVisitor {
    type Value = SerializableHandle;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("A Handle")
    }

    fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
        where
            E: de::Error,
    {
        Ok(SerializableHandle(Handle::from_integer(value)))
    }
}

impl<'de> Deserialize<'de> for SerializableHandle {
    fn deserialize<D>(deserializer: D) -> Result<SerializableHandle, D::Error>
        where
            D: Deserializer<'de>,
    {
        deserializer.deserialize_u64(SerializableHandleVisitor)
    }
}


// Implement both Deref and DerefMut for transparent usage
impl Deref for SerializableHandle {
    type Target = Handle;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SerializableHandle {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}