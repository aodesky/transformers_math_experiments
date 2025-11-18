/*
 * Compute log(conductors) in parallel using PARI/GP
 *
 * Reads data/ainvs_of_big_starter_dataset.txt and computes
 * the log of the conductor for each elliptic curve using parallel processing.
 *
 * Input format: rational,a1,a2,a3,a4,a6
 * Output format: rational,log_conductor
 */

default(threadsizemax, 100000000)

/* Function to compute log(conductor) from ainvs - must be self-contained for parapply */
compute_conductor(ainvs_data) =
{
    my(rational, ainvs, E, red, conductor, log_conductor);

    /* ainvs_data is a vector [rational, a1, a2, a3, a4, a6] */
    rational = ainvs_data[1];
    ainvs = [ainvs_data[2], ainvs_data[3], ainvs_data[4], ainvs_data[5], ainvs_data[6]];

    E = ellinit(ainvs);
    red = ellglobalred(E);
    conductor = red[1];
    log_conductor = log(conductor);

    /* Return [rational, log_conductor] as a vector */
    [rational, log_conductor]
}

print("Reading data/ainvs_of_big_starter_dataset.txt...");

/* Read the entire file */
lines = readstr("../data/ainvs_of_big_starter_dataset.txt");
n = length(lines);

print("Found ", n, " curves");
print("Parsing input data...");

/* Parse each line into [rational, a1, a2, a3, a4, a6] */
parsed_data = vector(n);
for(j=1, n, my(parts = strsplit(lines[j], ",")); parsed_data[j] = [parts[1], eval(parts[2]), eval(parts[3]), eval(parts[4]), eval(parts[5]), eval(parts[6])]);

print("Parsed ", n, " lines");
print("Computing log(conductors) in parallel (using ", default(nbthreads), " threads)...");
print("");

/* Apply compute_conductor in parallel to all curves */
results = parapply(compute_conductor, parsed_data);

print("");
print("Writing results to data/log_conductors_of_big_starter_dataset.txt...");

/* Write results to file */
write1("../data/log_conductors_of_big_starter_dataset.txt", "");  /* Clear file */

for(j=1, n, write("../data/log_conductors_of_big_starter_dataset.txt", results[j][1], ",", results[j][2]));

print("");
print("Done! Processed ", n, " curves.");

/* Optional: show some statistics */
log_conductors = vector(n, k, results[k][2]);
print("");
print("Statistics:");
print("  Minimum log(conductor): ", vecmin(log_conductors));
print("  Maximum log(conductor): ", vecmax(log_conductors));

quit;
