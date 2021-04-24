use handlegraph::handle::Handle;
use serde::{Serialize, Deserialize};
use std::ops::{Deref, DerefMut};

/// Create a wrapper for Handle so that it can be serialized
//#[derive(Serialize, Deserialize)]
struct SerializableHandle(Handle);

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