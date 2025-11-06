/*
 * SpaceNav.cpp
 *
 * An interplanetary mission planner that uses NASA's CSPICE toolkit
 * for ephemeris data and a local implementation of a Lambert solver
 * for trajectory calculation.
 *
 * Now includes two modes:
 * 1. Specific Date Analysis: User provides all dates.
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
 * @brief Mode 1: Analyzes a single, user-defined trajectory.
 */
void runSpecificDateMode(SpiceDouble sun_gm)
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

    double dv_capability = tsiolkovsky(dry_mass, fuel_mass, isp);

    try
    {
        TransferAnalysis transfer = calculate_transfer_dv(sun_gm, origin_body, departure_et, dest_body, arrival_et);

        // --- Print Results ---
        std::cout << "\n--- Mission Analysis ---" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        
        std::cout << "Time of Flight:       " << transfer.time_of_flight_days << " days" << std::endl;
        std::cout << "Required dV (Depart): " << transfer.dv_departure << " km/s" << std::endl;
        std::cout << "Required dV (Arrive): " << transfer.dv_arrival << " km/s" << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "Total Required dV:    " << transfer.dv_total << " km/s" << std::endl;
        std::cout << "Spacecraft Max dV:    " << dv_capability << " km/s" << std::endl;

        std::cout << "\n--- RESULT ---" << std::endl;
        if (dv_capability >= transfer.dv_total)
        {
            std::cout << "SUCCESS: Mission is feasible." << std::endl;
            std::cout << "Delta-V Margin: " << (dv_capability - transfer.dv_total) << " km/s" << std::endl;
        }
        else
        {
            std::cout << "FAILURE: Insufficient Delta-V." << std::endl;
            std::cout << "Delta-V Shortfall: " << (transfer.dv_total - dv_capability) << " km/s" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nMission Calculation Error: " << e.what() << std::endl;
        std::cerr << "This mission profile (dates/planets) is likely impossible." << std::endl;
    }
}

/**
 * @brief Mode 2: Searches for the minimum-energy transfer within a window.
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

    std::cout << "Enter Departure Date (e.g., 2035-04-25 12:00): ";
    std::getline(std::cin, dep_date);

    // --- Calculations ---
    SpiceDouble departure_et;
    str2et_c(dep_date.c_str(), &departure_et);
    checkSpiceError("str2et_c (Parsing departure date)");

    double dv_capability = tsiolkovsky(dry_mass, fuel_mass, isp);

    std::cout << "\n--- Searching for Optimal Transfer ---" << std::endl;
    std::cout << "Spacecraft Max dV: " << std::fixed << std::setprecision(3) << dv_capability << " km/s" << std::endl;
    std::cout << "Searching arrival window (Depart + 120 to + 365 days)..." << std::endl;

    // Define the search window
    const int search_start_days = 120; // Start searching 4 months out
    const int search_end_days   = 365; // Stop searching 1 year out
    
    double min_dv_required = std::numeric_limits<double>::max();
    SpiceDouble best_arrival_et = 0;
    TransferAnalysis best_transfer;

    for (int days = search_start_days; days <= search_end_days; ++days)
    {
        SpiceDouble arrival_et = departure_et + (days * DAY_IN_SECONDS);
        
        try
        {
            TransferAnalysis current_transfer = calculate_transfer_dv(sun_gm, origin_body, departure_et, dest_body, arrival_et);
            
            if (current_transfer.dv_total < min_dv_required)
            {
                min_dv_required = current_transfer.dv_total;
                best_arrival_et = arrival_et;
                best_transfer = current_transfer;
            }
        }
        catch (const std::exception& e)
        {
            // This transfer (e.g., day 121) was impossible.
            // We just ignore it and continue the loop.
        }
    }

    // --- Print Results ---
    if (best_arrival_et == 0)
    {
        std::cerr << "\n--- RESULT ---" << std::endl;
        std::cerr << "FAILURE: No possible transfer found in the 1-year search window." << std::endl;
        return;
    }

    // Convert the best arrival ET back to a string
    SpiceChar best_arr_date_str[128];
    et2utc_c(best_arrival_et, "ISOC", 3, 128, best_arr_date_str);
    checkSpiceError("et2utc_c (Converting best arrival time)");

    std::cout << "\n--- Optimal Transfer Found ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
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
    std::cout << "2. Find optimal transfer (Depart date + Fuel)" << std::endl;
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
