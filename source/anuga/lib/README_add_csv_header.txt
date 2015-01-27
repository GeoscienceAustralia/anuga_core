The function defined here allows programmatic addition of a CSV
header line to an existing CSV file (without header).  There
is a check that the header line to be added has the same number
of fields as the existing CSV file.

There is a 'be_green' parameter designed to change the copy algorithm
so that less memory is used, but the function will be slower.
This option is not yet implemented.
