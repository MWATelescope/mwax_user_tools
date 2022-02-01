This is a tool to undo the data block ordering in one or more MWA subfiles to generate a fully
(re) interleaved data stream, which can then be piped or passed on to subsequent
tools.

The tool optionally outputs the metadata to different files if these need to be collected

Example usage: For a directory with many subfiles in it, e.g.:
1321838416_1321838416_125.sub
1321838416_1321838424_125.sub
1321838416_1321838432_125.sub
1321838416_1321838440_125.sub
1321838416_1321838448_125.sub

If these files are named sensibly so that they are time-ordered by the shell, then they can be processed simply as:

cat *.sub | reinterleave_subfile | my_downstream_tool ....

to output the headers (all to one file)

cat *.sub | reinterleave_subfile -h header_dump_file.txt | my_downstream_tool

if you really want to write the samples to a file:

cat *.sub | reinterleave_subfile -o my_big_file_of_samples.dat

