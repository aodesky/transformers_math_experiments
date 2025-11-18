/* Compute log(conductors) in parallel - Interactive GP version */

/* Function to compute log(conductor) - must be self-contained for parapply */
compute_conductor(ainvs_data) = {
    my(rational = ainvs_data[1], ainvs = [ainvs_data[2], ainvs_data[3], ainvs_data[4], ainvs_data[5], ainvs_data[6]], E = ellinit(ainvs), red = ellglobalred(E), conductor = red[1], log_conductor = log(conductor));
    [rational, log_conductor]
}

/* Read and parse data */
print("Reading data/ainvs_of_big_starter_dataset.txt...");
lines = readstr("../data/ainvs_of_big_starter_dataset.txt");
n = length(lines);
print("Found ", n, " curves");

print("Parsing input data...");
parsed_data = vector(n);
for(j=1, n, my(parts = strsplit(lines[j], ",")); parsed_data[j] = [parts[1], eval(parts[2]), eval(parts[3]), eval(parts[4]), eval(parts[5]), eval(parts[6])]);
print("Parsed ", n, " lines");

/* Compute in parallel */
print("Computing log(conductors) in parallel (using ", default(nbthreads), " threads)...");
results = parapply(compute_conductor, parsed_data);
print("Done computing!");

/* Write results */
print("Writing results to data/log_conductors_of_big_starter_dataset.txt...");
write1("../data/log_conductors_of_big_starter_dataset.txt", "");
for(j=1, n, write("../data/log_conductors_of_big_starter_dataset.txt", results[j][1], ",", results[j][2]));
print("Done writing!");

/* Statistics */
log_conductors = vector(n, k, results[k][2]);
print("Statistics:");
print("  Minimum log(conductor): ", vecmin(log_conductors));
print("  Maximum log(conductor): ", vecmax(log_conductors));
print("Complete!");
