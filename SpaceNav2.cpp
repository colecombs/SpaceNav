/*
 * SpaceNav.cpp
 *
 * An interplanetary mission planner that uses NASA's CSPICE toolkit
 * for ephemeris data and a local implementation of a Lambert solver
 * for trajectory calculation.
 *
 * Now includes two modes:
 * 1. Specific Date Analysis: User provides all dates, fuel is calculated.
 * 2. Transfer Search: User provides a departure date, and the
 * program searches for the minimum-energy transfer.
 */

// Define _USE_MATH_DEFINES to get M_PI from <cmath>
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip> // For std::setw and std::setprecision
#include <limits>  // For std::numeric_limits
#include <stdexcept> // For std::runtime_error

// Include the SPICE user toolkit
extern "C"
{
#include "SpiceUsr.h"
}

// Include our local header-only Lambert solver
#include "lambert_solver.h"

// --- Global Constants ---
const double DAY_IN_SECONDS = 86400.0;
// Standard gravity in m/s^2. We divide by 1000 to get km/s^2.
const double G0_KMS = 9.80665 / 1000.0;


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
double calculate_dv(double dry_mass, double fuel_mass, double isp)
{
    double total_mass = dry_mass + fuel_mass;
    
    // Tsiolkovsky rocket equation: dV = Isp * g0 * ln(m_total / m_dry)
    double delta_v = isp * G0_KMS * std::log(total_mass / dry_mass);
    
    return delta_v;
}

/**
 * @brief Calculates the required fuel mass for a given Delta-V.
 * This is the inverted Tsiolkovsky rocket equation.
 * @param dry_mass The mass of the spacecraft without fuel (kg).
 * @param delta_v The required change in velocity (km/s).
 * @param isp The specific impulse of the engine (seconds).
 * @return The required propellant mass (kg).
 */
double calculate_fuel(double dry_mass, double delta_v, double isp)
{
    // m_total = m_dry * exp(dV / (Isp * g0))
    double total_mass = dry_mass * std::exp(delta_v / (isp * G0_KMS));
    
    // fuel_mass = m_total - m_dry
    double fuel_mass = total_mass - dry_mass;
    
    return fuel_mass;
}


/**
 * @brief A helper struct to hold the results of a transfer calculation.
 */
struct TransferAnalysis
{
    double dv_departure;
    double dv_arrival;
    double dv_total;
    double time_of_flight_days;
};

/**
 * @brief Encapsulates the core SPICE and Lambert logic.
 * Throws a runtime_error if the transfer is impossible.
 *
 * @param sun_gm Gravitational parameter of the Sun.
 * @param origin_body Name of the origin body.
 * @param departure_et Ephemeris time of departure.
 *TA
 * @param dest_body Name of the destination body.
 * @param arrival_et Ephemeris time of arrival.
 * @return A TransferAnalysis struct with the calculated dV values.
 */
TransferAnalysis calculate_transfer_dv(SpiceDouble sun_gm, const std::string& origin_body, SpiceDouble departure_et, const std::string& dest_body, SpiceDouble arrival_et)
{
    // --- 6. Get Planet States from SPICE ---
    SpiceDouble origin_state[6], dest_state[6];
    SpiceDouble light_time;

    spkezr_c(origin_body.c_str(), departure_et, "J2000", "NONE", "SUN", origin_state, &light_time);
    checkSpiceError("spkezr_c for origin");
    
    spkezr_c(dest_body.c_str(), arrival_et, "J2000", "NONE", "SUN", dest_state, &light_time);
    checkSpiceError("spkezr_c for destination");

    // --- 7. Solve Lambert's Problem ---
    std::vector<double> r1 = {origin_state[0], origin_state[1], origin_state[2]};
    std::vector<double> r2 = {dest_state[0], dest_state[1], dest_state[2]};
    std::vector<double> v1_transfer(3);
    std::vector<double> v2_transfer(3);
    
    double time_of_flight_seconds = arrival_et - departure_et;
    if (time_of_flight_seconds <= 0)
    {
        throw std::runtime_error("Arrival date is before departure date.");
    }

    // Call our local Lambert solver (throws error on failure)
    LambertSolver::solve_lambert(sun_gm, r1.data(), r2.data(), time_of_flight_seconds, true, v1_transfer.data(), v2_transfer.data());

    // --- 8. Calculate Required Delta-V ---
    TransferAnalysis result;
    result.time_of_flight_days = time_of_flight_seconds / DAY_IN_SECONDS;

    double v_inf_dep[3];
    vsub_c(v1_transfer.data(), &origin_state[3], v_inf_dep);
    result.dv_departure = vnorm_c(v_inf_dep);

    double v_inf_arr[3];
    vsub_c(v2_transfer.data(), &dest_state[3], v_inf_arr);
    result.dv_arrival = vnorm_c(v_inf_arr);

    result.dv_total = result.dv_departure + result.dv_arrival;
    
    return result;
}


/**
 * @brief Mode 1: Analyzes a single, user-defined trajectory and calculates required fuel.
 */
void runSpecificDateMode(SpiceDouble sun_gm)
{
    // --- Get User Input ---
    double dry_mass, isp;
    std::cout << "\n--- Spacecraft Parameters ---" << std::endl;
    std::cout << "Enter Dry Mass (kg): ";
    std::cin >> dry_mass;
    
    std::cout << "Enter Engine ISP (seconds): ";
    std::cin >> isp;

    std::string origin_body, dest_body, dep_date, arr_date;
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

    // --- Calculations ---
    SpiceDouble departure_et, arrival_et;
    str2et_c(dep_date.c_str(), &departure_et);
    checkSpiceError("str2et_c (Parsing departure date)");
    
    str2et_c(arr_date.c_str(), &arrival_et);
    checkSpiceError("str2et_c (Parsing arrival date)");

    try
    {
        TransferAnalysis transfer = calculate_transfer_dv(sun_gm, origin_body, departure_et, dest_body, arrival_et);
        double required_fuel = calculate_fuel(dry_mass, transfer.dv_total, isp);

        // --- Print Results ---
        std::cout << "\n--- Mission Analysis ---" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        
        std::cout << "Time of Flight:       " << transfer.time_of_flight_days << " days" << std::endl;
        std::cout << "Required dV (Depart): " << transfer.dv_departure << " km/s" << std::endl;
        std::cout << "Required dV (Arrive): " << transfer.dv_arrival << " km/s" << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "Total Required dV:    " << transfer.dv_total << " km/s" << std::endl;
        
        std::cout << "\n--- RESULT ---" << std::endl;
        std::cout << "Required Fuel Mass:   " << required_fuel << " kg" << std::endl;
        std::cout << "Total Spacecraft Mass:" << (dry_mass + required_fuel) << " kg" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nMission Calculation Error: " << e.what() << std::endl;
        std::cerr << "This mission profile (dates/planets) is likely impossible." << std::endl;
    }
}

/**
 * @brief Mode 2: Searches for the minimum-energy transfer window (2D "Porkchop Plot" search).
 */
void runSearchMode(SpiceDouble sun_gm)
{
    // --- Get User Input ---
    double dry_mass, fuel_mass, isp;
    std::cout << "\n--- Spacecraft Parameters ---" << std::endl;
    std::cout << "Enter Dry Mass (kg): ";
    std::cin >> dry_mass;
    
    std::cout << "Enter Fuel Mass (kg): ";
    std::cin >> fuel_mass;
    
    std::cout << "Enter Engine ISP (seconds): ";
    std::cin >> isp;

    std::string origin_body, dest_body, dep_date;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::cout << "\n--- Mission Parameters ---" << std::endl;
    std::cout << "Enter Origin Body (e.g., EARTH): ";
    std::getline(std::cin, origin_body);

    std::cout << "Enter Destination Body (e.g., MARS BARYCENTER): ";
    std::getline(std::cin, dest_body);

    std::cout << "Enter *Start* of Departure Window (e.g., 2035-01-01): ";
    std::getline(std::cin, dep_date);

    // --- Calculations ---
    SpiceDouble departure_start_et;
    str2et_c(dep_date.c_str(), &departure_start_et);
    checkSpiceError("str2et_c (Parsing departure start date)");

    double dv_capability = calculate_dv(dry_mass, fuel_mass, isp);

    std::cout << "\n--- Searching for Optimal Transfer Window ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Spacecraft Max dV: " << dv_capability << " km/s" << std::endl;

    // --- Define the 2D Search Space ---
    // We will search a 1-year departure window
    const int departure_window_days = 365; 
    // For each departure, we will check a flight time window
    const int tof_min_days = 120; // Min time of flight
    const int tof_max_days = 400; // Max time of flight
    // Search step size. 5 days is faster. 1 day is more precise.
    const int step_days = 5; 

    std::cout << "Searching Departure Window: " << dep_date << " to ~1 year after." << std::endl;
    std::cout << "Searching Flight Times: " << tof_min_days << " to " << tof_max_days << " days" << std::endl;
    std::cout << "This may take a moment..." << std::endl;

    double min_dv_required = std::numeric_limits<double>::max();
    SpiceDouble best_departure_et = 0;
    SpiceDouble best_arrival_et = 0;
    TransferAnalysis best_transfer;

    // --- Begin 2D Search (Nested Loop) ---
    for (int dep_day = 0; dep_day <= departure_window_days; dep_day += step_days)
    {
        SpiceDouble current_departure_et = departure_start_et + (dep_day * DAY_IN_SECONDS);

        for (int tof_day = tof_min_days; tof_day <= tof_max_days; tof_day += step_days)
        {
            SpiceDouble current_arrival_et = current_departure_et + (tof_day * DAY_IN_SECONDS);
            
            try
            {
                TransferAnalysis current_transfer = calculate_transfer_dv(sun_gm, origin_body, current_departure_et, dest_body, current_arrival_et);
                
                if (current_transfer.dv_total < min_dv_required)
                {
                    // This is our new best trajectory
                    min_dv_required = current_transfer.dv_total;
                    best_departure_et = current_departure_et;
                    best_arrival_et = current_arrival_et;
                    best_transfer = current_transfer;
                }
            }
            catch (const std::exception& e)
            {
                // This combination of dates is impossible.
                // We just ignore it and continue the loop.
            }
        } // end Time-of-Flight loop
    } // end Departure Date loop

    // --- Print Results ---
    if (best_arrival_et == 0)
    {
        std::cerr << "\n--- RESULT ---" << std::endl;
        std::cerr << "FAILURE: No possible transfer found in the search window." << std::endl;
        return;
    }

    // Convert the best dates back to strings
    SpiceChar best_dep_date_str[128];
    SpiceChar best_arr_date_str[128];
    et2utc_c(best_departure_et, "ISOC", 3, 128, best_dep_date_str);
    checkSpiceError("et2utc_c (Converting best departure time)");
    et2utc_c(best_arrival_et, "ISOC", 3, 128, best_arr_date_str);
    checkSpiceError("et2utc_c (Converting best arrival time)");

    std::cout << "\n--- Optimal Transfer Window Found ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    std::cout << "Best Departure Date:  " << best_dep_date_str << " UTC" << std::endl;
    std::cout << "Best Arrival Date:    " << best_arr_date_str << " UTC" << std::endl;
    std::cout << "Time of Flight:       " << best_transfer.time_of_flight_days << " days" << std::endl;
    std::cout << "Required dV (Depart): " << best_transfer.dv_departure << " km/s" << std::endl;
    std::cout << "Required dV (Arrive): " << best_transfer.dv_arrival << " km/s" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Total Required dV:    " << best_transfer.dv_total << " km/s" << std::endl;
    std::cout << "Spacecraft Max dV:    " << dv_capability << " km/s" << std::endl;

    std::cout << "\n--- RESULT ---" << std::endl;
    if (dv_capability >= best_transfer.dv_total)
    {
        std::cout << "SUCCESS: Mission is feasible with this optimal trajectory." << std::endl;
        std::cout << "Delta-V Margin: " << (dv_capability - best_transfer.dv_total) << " km/s" << std::endl;
    }
    else
    {
        std::cout << "FAILURE: Insufficient Delta-V for the *most efficient* transfer." << std::endl;
        std::cout << "Delta-V Shortfall: " << (best_transfer.dv_total - dv_capability) << " km/s" << std::endl;
    }
}


// --- Main Program ---
int main()
{
    // --- 1. Load SPICE Kernels ---
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
    
    // --- 3. Mode Selection ---
    int mode = 0;
    std::cout << "\n--- Select Mode ---" << std::endl;
    std::cout << "1. Analyze a specific trajectory (Depart + Arrive dates)" << std::endl;
    std::cout << "2. Find optimal transfer window (Depart date + Fuel)" << std::endl;
    std::cout << "Enter mode (1 or 2): ";
    std::cin >> mode;

    switch (mode)
    {
        case 1:
            runSpecificDateMode(sun_gm);
            break;
        case 2:
            runSearchMode(sun_gm);
            break;
        default:
            std::cerr << "Invalid mode selected. Exiting." << std::endl;
            return 1;
    }
    
    return 0;
}
