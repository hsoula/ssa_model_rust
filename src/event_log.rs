use std::collections::HashMap;
use crate::reactions_maps::get_reaction_name;
pub struct Event {
    time: f64,
    pub even_type: u32
}


pub struct EventLog {
    events: Vec<Event>,
}
impl Event {
    fn new(time: f64, eventtype: u32) -> Event {
        Event {time:time, even_type:eventtype}
    }
}

impl EventLog {
    pub fn new() -> EventLog {
        EventLog {events:Vec::new()}
    }

    pub fn get_number_of_event(&self, event_type: u32) -> u32{
        self.events.iter()
            .filter(|e| e.even_type == event_type)
            .collect::<Vec<&Event>>().len() as u32
    }

    pub fn event_stats(&self) -> HashMap<String,i32> {
        let event_states = self.events.iter()
            .map(|event| get_reaction_name(event.even_type).unwrap())
            .fold(HashMap::new(), |mut acc, key| {
                *acc.entry(key).or_insert(0) += 1;
                acc
            });
        event_states
    }
    pub fn add_event(&mut self, time: f64, event_type: u32) {
        let event = Event::new(time, event_type);
        self.events.push(event);
    }
}