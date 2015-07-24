#!/bin/bash

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


