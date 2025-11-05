/*
 * SpaceNav.cpp
 *
 * An interplanetary mission planner that uses NASA's CSPICE toolkit
 * for ephemeris data and a local implementation of a Lambert solver
 * for trajectory calculation.
 */

// Define _USE_MATH_DEFINES to get M_PI from <cmath>
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip> // For std::setw and std::setprecision
#include <limits>  // For std::numeric_limits

// Include the SPICE user toolkit
extern "C"
{
#include "SpiceUsr.h"
}

// Include our local header-only Lambert solver
#include "lambert_solver.h"

/**
 * @brief Utility function to check for SPICE errors after a call.
 * If an error is found, it prints the long message and exits.
 */
void checkSpiceError(const char *operation)
{
    if (failed_c())
    {
        std::cerr << "SPICE Error detected during: " << operation << std::endl;
        
        // Get and print the long error message
        SpiceChar longMessage[1024];
        getmsg_c("LONG", 1024, longMessage);
        std::cerr << longMessage << std::endl;

        // Reset the error status and exit
        reset_c();
        exit(1);
    }
}

/**
 * @brief Calculates the maximum Delta-V of a spacecraft using the
 * Tsiolkovsky rocket equation.
 * @param dry_mass The mass of the spacecraft without fuel (kg).
 * @param fuel_mass The mass of the propellant (kg).
 * @param isp The specific impulse of the engine (seconds).
 * @return The maximum Delta-V (km/s).
 */
double tsiolkovsky(double dry_mass, double fuel_mass, double isp)
{
    // Standard gravity in m/s^2. We divide by 1000 to get km/s^2.
    const double g0 = 9.80665 / 1000.0;
    
    double total_mass = dry_mass + fuel_mass;
    
    // Tsiolkovsky rocket equation: dV = Isp * g0 * ln(m_total / m_dry)
    double delta_v = isp * g0 * std::log(total_mass / dry_mass);
    
    return delta_v;
}


// --- Main Program ---
int main()
{
    // --- 1. Load SPICE Kernels ---
    // We use a "meta-kernel" (a .mk file) to list all the
    // data files we want to load.
    furnsh_c("kernels.mk");
    checkSpiceError("furnsh_c (Loading meta-kernel)");

    // --- 2. Get Toolkit Version and Sun's Gravitational Parameter ---
    const SpiceChar *tkVersion = tkvrsn_c("TOOLKIT");
    checkSpiceError("tkvrsn_c (Getting toolkit version)");
    
    std::cout << "--- Interplanetary Mission Planner ---" << std::endl;
    std::cout << "Using SPICE Toolkit Version: " << tkVersion << std::endl;

    SpiceDouble sun_gm;
    SpiceInt n_gm;
    bodvcd_c(10, "GM", 1, &n_gm, &sun_gm);
    checkSpiceError("bodvcd_c (Getting Sun's GM)");
    
    // --- 3. Get User Input for Spacecraft ---
    double dry_mass, fuel_mass, isp;
    std::cout << "\n--- Spacecraft Parameters ---" << std::endl;
    std::cout << "Enter Dry Mass (kg): ";
    std::cin >> dry_mass;
    
    std::cout << "Enter Fuel Mass (kg): ";
    std::cin >> fuel_mass;
    
    std::cout << "Enter Engine ISP (seconds): ";
    std::cin >> isp;

    // --- 4. Get User Input for Mission ---
    std::string origin_body, dest_body, dep_date, arr_date;

    // We must "ignore" the newline character left in the input buffer
    // after the last `std::cin >> isp;`. This is critical for 
    // `std::getline` to work correctly.
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::cout << "\n--- Mission Parameters ---" << std::endl;
    std::cout << "Enter Origin Body (e.g., EARTH): ";
    std::getline(std::cin, origin_body);

    std::cout << "Enter Destination Body (e.g., MARS BARYCENTER): ";
    std::getline(std::cin, dest_body);

    std::cout << "Enter Departure Date (e.g., 2035-04-25 12:00): ";
    std::getline(std::cin, dep_date);

    std::cout << "Enter Arrival Date (e.g., 2036-01-08 09:00): ";
    std::getline(std::cin, arr_date);

    // --- 5. Convert Times to Ephemeris Time (ET) ---
    SpiceDouble departure_et, arrival_et;
    str2et_c(dep_date.c_str(), &departure_et);
    checkSpiceError("str2et_c (Parsing departure date)");
    
    str2et_c(arr_date.c_str(), &arrival_et);
    checkSpiceError("str2et_c (Parsing arrival date)");

    double time_of_flight_seconds = arrival_et - departure_et;
    if (time_of_flight_seconds <= 0)
    {
        std::cerr << "Error: Arrival date must be after departure date." << std::endl;
        return 1;
    }

    // --- 6. Get Planet States from SPICE ---
    // We get the state (position + velocity) of the planets relative
    // to the Sun ("SSB" or 10) in the "J2000" reference frame.
    // "NONE" specifies no light-time correction.
    SpiceDouble origin_state[6], dest_state[6];
    SpiceDouble light_time; // Unused, but required by spkezr_c

    spkezr_c(origin_body.c_str(), departure_et, "J2000", "NONE", "SUN", origin_state, &light_time);
    checkSpiceError("spkezr_c for origin");
    
    spkezr_c(dest_body.c_str(), arrival_et, "J2000", "NONE", "SUN", dest_state, &light_time);
    checkSpiceError("spkezr_c for destination");

    // --- 7. Solve Lambert's Problem ---
    // This finds the velocity vectors at the start and end of a transfer orbit
    // connecting two position vectors in a given time.
    
    // Extract position vectors (first 3 elements of the state vector)
    std::vector<double> r1 = {origin_state[0], origin_state[1], origin_state[2]};
    std::vector<double> r2 = {dest_state[0], dest_state[1], dest_state[2]};
    
    // Output velocity vectors
    std::vector<double> v1_transfer(3);
    std::vector<double> v2_transfer(3);

    // Call our local Lambert solver
    // We assume a prograde, 0-revolution transfer
    bool prograde = true;
    
    try 
    {
        // Pass sun_gm first, and use .data() to pass the C-style array
        LambertSolver::solve_lambert(sun_gm, r1.data(), r2.data(), time_of_flight_seconds, prograde, v1_transfer.data(), v2_transfer.data());
    } 
    catch (const std::exception& e)
    {
        std::cerr << "\nLambert Solver Error: " << e.what() << std::endl;
        std::cerr << "This mission profile (dates/planets) is likely impossible. Try different dates." << std::endl;
        return 1;
    }

    // --- 8. Calculate Required Delta-V ---
    // The required velocity change is the difference between the spacecraft's
    // orbital velocity and the required transfer orbit velocity.

    // v_infinity_dep = v_transfer - v_planet_at_departure
    double v_inf_dep[3];
    vsub_c(v1_transfer.data(), &origin_state[3], v_inf_dep);
    double dv_departure = vnorm_c(v_inf_dep);

    // v_infinity_arr = v_transfer - v_planet_at_arrival
    double v_inf_arr[3];
    vsub_c(v2_transfer.data(), &dest_state[3], v_inf_arr);
    double dv_arrival = vnorm_c(v_inf_arr);

    double dv_total_required = dv_departure + dv_arrival;

    // --- 9. Calculate Spacecraft Capability & Print Results ---
    double dv_capability = tsiolkovsky(dry_mass, fuel_mass, isp);

    std::cout << "\n--- Mission Analysis ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3); // Format output to 3 decimal places
    
    std::cout << "Time of Flight:       " << time_of_flight_seconds / 86400.0 << " days" << std::endl;
    std::cout << "Required dV (Depart): " << dv_departure << " km/s" << std::endl;
    std::cout << "Required dV (Arrive): " << dv_arrival << " km/s" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Total Required dV:    " << dv_total_required << " km/s" << std::endl;
    std::cout << "Spacecraft Max dV:    " << dv_capability << " km/s" << std::endl;

    std::cout << "\n--- RESULT ---" << std::endl;
    if (dv_capability >= dv_total_required)
    {
        std::cout << "SUCCESS: Mission is feasible." << std::endl;
        std::cout << "Delta-V Margin: " << (dv_capability - dv_total_required) << " km/s" << std::endl;
    }
    else
    {
        std::cout << "FAILURE: Insufficient Delta-V." << std::endl;
        std::cout << "Delta-V Shortfall: " << (dv_total_required - dv_capability) << " km/s" << std::endl;
    }
    
    return 0;
}

