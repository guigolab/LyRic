#!/bin/bash

sed 's/.*\t//g'| sed 's/; /\n/g'| sed 's/=/\t/g'
