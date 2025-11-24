pub struct Link {
    pub id : u32,
    pub links : [u32; 2],
    pub destroy : bool
}

impl Link {
    pub fn new(id:u32,a:u32, b:u32) -> Link {
        Link{id: id, links: [a, b], destroy:false}
    }
}