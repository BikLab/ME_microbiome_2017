#!/bin/bash

grep ">" $1 | cut -f 1 -d " " > $2
