#!/usr/bin/env python3

OUTPUT_INPUT_FILE = "set_def.in"
OUTPUT_OBS_FILE = "set_def.out.JLA_MWR12"

U, V, P, T = 1, 2, 4, 5

# Test configuration: 1 longitude, 9 latitudes, 1 model level
longitudes = [3 + 18 * i for i in range(20)]
latitudes = [-72, -54, -36, -18, 0, 18, 36, 54, 72]
levels = range(1, 6)

NOBS = (
    3 * len(levels) * len(latitudes) * len(longitudes)
    + 1 * len(latitudes) * len(longitudes)
)

def write_3d_obs(f, obs_type, variance):
    """Write observations at all horizontal locations and model levels."""
    for level in levels:
        for lat in latitudes:
            for lon in longitudes:
                f.write(
                    f"0\n"
                    f"{obs_type}\n"
                    f"1\n"
                    f"{level}\n"
                    f"{lon}\n"
                    f"{lat}\n"
                    f"0 0\n"
                    f"{variance}\n"
                )


def write_surface_pressure_obs(f, variance):
    """Write surface-pressure observations at all horizontal locations."""
    for lat in latitudes:
        for lon in longitudes:
            f.write(
                f"0\n"
                f"{P}\n"
                f"-1\n"
                f"-100\n"
                f"{lon}\n"
                f"{lat}\n"
                f"0 0\n"
                f"{variance}\n"
            )


with open(OUTPUT_INPUT_FILE, "w") as f:
    # Header
    f.write(f"{NOBS}\n")
    f.write("0\n")
    f.write("0\n")

    write_3d_obs(f, U, 9)
    write_3d_obs(f, V, 9)
    write_surface_pressure_obs(f, 40000)
    write_3d_obs(f, T, 9)

    # Once create_obs_sequence has read NOBS observations,
    # it asks directly for the output filename.
    f.write(f"{OUTPUT_OBS_FILE}\n")

print(f"Created {OUTPUT_INPUT_FILE} with {NOBS} observations.")
