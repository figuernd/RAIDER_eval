#!/bin/bash

# create_pals.sh: given a string $rep, creates all seeds that are the result of
# continuously concatenating $rep with itself. Puts seeds of length greater than
# specified by $min but less than specified by $max into file $fname
# Arguments: rep min max fname
# by Carly Schaeffer


#rep is the string to be repeated
str=${rep}
until [ "${#str}" -ge "${min}" ]; do
    str=${str}${rep}
done

until [ ${#str} -ge ${max} ]; do
    strcur=${str}1
    echo "${strcur}" >> $fname
    str=${str}${rep}
done
