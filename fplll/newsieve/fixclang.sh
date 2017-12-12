#!/bin/bash
clang-format -style=file $1 > CLANGTMP.txt
meld $1 CLANGTMP.txt
