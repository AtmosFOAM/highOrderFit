#!/bin/bash
set -e
wmake -j -a src
wmake -j test
