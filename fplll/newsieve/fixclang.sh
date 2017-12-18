#!/bin/bash
clang-format -style=file $1 > CLANGTMP.txt
diff $1 CLANGTMP.txt
