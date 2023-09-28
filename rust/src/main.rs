use anyhow::{Context, Result};
use chrono::{Datelike, TimeZone, Timelike};
use chrono_tz::UTC;
use std::fs;

fn main() -> Result<()> {
    // Load the OMM data from the JSON file
    let data = fs::read_to_string("space-track-omm.json").expect("Unable to read file");
    let satellites: Vec<sgp4::Elements> =
        serde_json::from_str(&data).expect("JSON was not well-formatted");
    println!("Loaded {} satellites", satellites.len());

    // Given coordinates
    let lat = 34.56;
    let lon = -118.76;

    // Get the current time
    let now = chrono::Utc::now();
    // Iterate over the satellites, propagate their orbits, and find the closest one
    let predictions = satellites
        .iter()
        .map(|sat| {
            let constants = sgp4::Constants::from_elements(sat).unwrap();
            let sat_ndt = sat.datetime;
            let sat_utc_dt = chrono::Utc.from_local_datetime(&sat_ndt).unwrap();
            let time_diff = now - sat_utc_dt;
            let epoch_minutes = (time_diff.num_seconds() as f64) / 60.0;
            println!("Epoch: {:?} {}", sat_utc_dt, epoch_minutes);
            let prediction = constants.propagate(epoch_minutes).unwrap();
            // The sgp4 docs say "The position and velocity are given in the
            // True Equator, Mean Equinox (TEME) of epoch reference frame" but
            // we need to convert to lat, lon, altitude.

        })
        .collect::<Vec<_>>();
    Ok(())
}

// Function to calculate the distance between two coordinates
fn haversine_distance(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let r = 6371.0; // Radius of the Earth in km
    let d_lat = (lat2 - lat1).to_radians();
    let d_lon = (lon2 - lon1).to_radians();
    let a = (d_lat / 2.0).sin().powi(2)
        + lat1.to_radians().cos() * lat2.to_radians().cos() * (d_lon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
    r * c
}
