# Find all lines in Python code which contain multiple old_div
for repeats in 1 2 3 4 5 ; do
   find ./ -name "*.py" -exec grep -EHn "(old_div.*){${repeats}}" {} \; >/tmp/old_div${repeats}
   echo $repeats "or more repeats of old_div per line: number of lines:" `wc -l /tmp/old_div${repeats}`
done

