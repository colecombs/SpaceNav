/*
 * HelloMars.cpp
 *
 * This is a "Hello, World!" program for the CSPICE toolkit.
 *
 * It demonstrates the basic structure of a SPICE program:
 * 1. Load necessary kernels (data files).
 * 2. Perform a calculation.
 * 3. Unload kernels (good practice, though not strictly needed here).
 *
 * This program calculates the position of Mars (499) relative to
 * Earth (399) in the J2000 reference frame.
 */

// We use C-style I/O (printf) for simple output
// and C++ streams (iostream) for error messages.
#include <iostream>
#include <stdio.h>

// The one and only C-style header file you need for CSPICE.
// It contains all function declarations.
extern "C" {
    #include "SpiceUsr.h"
}

// A simple helper function to check for SPICE errors
// after any 'cspice' call.
void checkSpiceError(const char* operation) {
    // If a SPICE error was signaled, 'failed_c()' will be true.
    if (failed_c()) {
        // We'll create a buffer to hold the long error message.
        // SPICE error messages can be quite long.
        SpiceChar longMessage[1024];
        
        // Get the long error message
        getmsg_c("LONG", 1024, longMessage);

        // Print it to std::cerr (the error stream)
        std::cerr << "SPICE Error detected during: " << operation << std::endl;
        std::cerr << longMessage << std::endl;

        // 'reset_c()' clears the error status
        // so we can continue (or, in a real app, exit).
        reset_c();

        // In a real application, you might want to exit:
        // exit(1);
    }
}


int main() {
    // --- 1. Load Kernels ---
    // A "metakernel" is a text file that lists all other
    // kernels we want to load. This is the best practice.
    //
    // We haven't created this file yet, but we will in the next step.
    const char* metakernel = "kernels.mk";
    furnsh_c(metakernel);
    checkSpiceError("furnsh_c (loading metakernel)");

    // --- 2. Define the Calculation ---
    SpiceDouble et;             // Ephemeris Time (the time of our query)
    SpiceDouble lightTime;      // One-way light time from target to observer
    SpiceDouble marsState[6];   // Output state vector (x,y,z, vx,vy,vz)
    
    // We'll ask for the state of Mars "now".
    // 'str2et_c' converts a time string to Ephemeris Time.
    //
    // EDIT: "now" is not a reliable string for str2et_c.
    // We will use a specific ISO 8601 time string, which is
    // guaranteed to work. Let's use the current time.
    const char* timeString = "2025-10-31T18:50:00";
    str2et_c(timeString, &et);
    checkSpiceError("str2et_c (converting time string)");

    // --- 3. Perform the Calculation ---
    // This is the core SPICE call.
    //
    // "MARS": The target body
    // "EARTH": The observing body
    // "J2000": The reference frame (our "map grid")
    // "NONE":  Aberration correction (we'll keep it simple)
    //
    // This function looks up the data from the kernels we loaded.
    //
    // FIX: The de440s.bsp kernel has the Earth-Moon Barycenter (3),
    // not "EARTH" (399) as a separate body. We'll use "3" as the observer.
    //
    // NEW FIX: The de421.bsp kernel *does* contain "EARTH" (399),
    // so we can go back to using "EARTH" as the observer.
    spkezr_c("MARS", et, "J2000", "NONE", "EARTH", marsState, &lightTime);
    checkSpiceError("spkezr_c (calculating Mars state)");

    // --- 4. Display Results ---
    // The first 3 elements of marsState are the x,y,z position
    // in kilometers.
    printf("\n");
    printf("Position of Mars relative to Earth at %s UTC:\n", timeString);
    printf("  X: %15.3f km\n", marsState[0]);
    printf("  Y: %15.3f km\n", marsState[1]);
    printf("  Z: %15.3f km\n", marsState[2]);
    printf("\n");
    printf("One-way light time: %15.3f seconds\n", lightTime);
    printf("\n");

    // --- 5. Unload Kernels ---
    // 'kclear_c' unloads all kernels. It's good hygiene.
    kclear_c();
    checkSpiceError("kclear_c (unloading kernels)");

    return 0;
}




