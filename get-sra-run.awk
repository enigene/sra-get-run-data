# Get SRA Run data v1.1, Nov 10 2013
# Download SRA sequences directly.
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# v1.0, Oct 12 2103 - Initial version
# v1.1, Nov 10 2103 - Added fastq load option. Kinda portable improvement.

# Usage: gawk -v sra=SRR385751 -v num=200 -v random=1 -v fastq=1 -f get-sra-run.awk

BEGIN {
  MAX_SPOT_PER_PAGE = 10000
  SERVER = "http://trace.ncbi.nlm.nih.gov/Traces/sra/"

  HIGH_CHAR = 93 # chr(93) == ~
  LOW_CHAR = 33  # chr(33) == !

# Check platform
# http://stackoverflow.com/questions/12025617/awk-find-out-if-running-on-windows
  is_windows = 0
#ENVIRON["OS"]
  if (index(tolower(ENV["OS"]), "windows") > 0) {
    is_windows = 1
  }

# Get the Run accession
  if (sra != "") {
    runAcc = sra
  } else {
    print "Error: Run accession not specified! Try -v sra=SRRnnnnnn"
    exit
  }

# Set record identifier
  if (fastq) {
    rID = "@"
  } else {
    rID = ">"
  }

# Testing curl
  if (system(curl) != 0) {
    print "Error: program \"cURL\" not found!"
    system("")
    exit
  }

# Testing sleep
  if (system("sleep 1") != 0) {
    print "Error: program \"sleep\" not found!"
    system("")
    exit
  }

# Testing server
  if (is_windows) {
    test = "curl -s -I -o NUL -w %%{http_code} --url " SERVER
  } else {
    test = "curl -s -I -o /dev/null -w %{http_code} --url " SERVER
  }
  test | getline response_code
  close(test)
  if (response_code != 200) {
    print "Error: server " SERVER " not OK (>.<) (response code " response_code ")"
    exit
  } else {
    print "Server " SERVER " OK (response code " response_code ")"
  }

# Try get one spot from this sra

  test_query = "curl -s \"" SERVER "sra.cgi?run_spot=" sra "&page=" 1 "&page_size=" 1 "\""

  while ((test_query | getline) > 0) {
#   Get the total number of spots
    if($0 ~ /^amount/) {
      if (match($0, /:([0-9]+),/, tmp)) totalSpots = tmp[1]
      if (totalSpots)
        totalPages = ceil(totalSpots / MAX_SPOT_PER_PAGE)
    }
  }
  close(test_query)
  system("sleep 2")

  if (totalSpots == 0) {
    print "Error: Run accession (" runAcc ") has zero entries"
    exit
  } else {
    print "Total spots: " totalSpots
    print "Total pages: " totalPages
    print "Max spots per page: " MAX_SPOT_PER_PAGE
  }

# Calculate pages and spots per page
  if (num > 0) {
    if (num > MAX_SPOT_PER_PAGE) {
      pageSize = gcd(num, MAX_SPOT_PER_PAGE)
      pages = num / pageSize
    }
    if (num < MAX_SPOT_PER_PAGE) {
      pageSize = num
      pages = 1
    }
    if (num == MAX_SPOT_PER_PAGE) {
      pageSize = MAX_SPOT_PER_PAGE
      pages = 1
    }
  } else {
    pageSize = MAX_SPOT_PER_PAGE
    pages = 1
    num = 1
  }

# Creating empty array
  delete(pagesA)

# Populate array with page numbers as index
  if (random) {
    srand()
    while (length(pagesA) < pages) {
      randomNumber = sprintf("%07i", rand() * totalPages)
      pagesA[randomNumber]
    }
  } else {
    for (i = 1; i <= pages; i++) pagesA[sprintf("%07i", i)]
  }

# Output filename
  if (fastq) {
    outFile = runAcc "-" num ".fastq"
  } else {
    outFile = runAcc "-" num ".fas"
  }

# Retrive Reads and sequences
  l = asorti(pagesA)
  for (p = 1; p <= l; p++) {
    currPage = pagesA[p]+0
    firstRID = ""

    query = "curl -s \"" SERVER "sra.cgi?run_spot=" sra "&page=" currPage "&page_size=" pageSize "\""

    while ((query | getline) > 0) {
#     Get the Spot index
      if ($0 ~ /^id/) {
        readID = 1
        if (match($0, /:(.+?),/, tmp)) runId = tmp[1]
        if (!firstRID) firstRID = runId
      }
#     Get Read sequence
      if ($0 !~ (/[^ATGCN",]/)) {
        if (match($0, /\"([ATGCN]+?)\"/, tmp)) read = tmp[1]
        if (runId && read) {
          pairA[readID] = rID "gnl|SRA|" runAcc "." runId "." readID "\n" read "\n"
          readID++
          j++
        }
        read = ""
      }
      if (fastq) {
        if ($0 ~ (/^qs:/)) { qs = 1; readID = 1 }
        if ($0 ~ /^\],/) { qs = 0 }
        if (qs) {
          if (match($0, /\[([[:digit:],]+?)\]/, tmp)) qual = tmp[1]
          if (qual) {
            qlen = split(qual, qualS, ",")
            qpos = 0
            while (qpos < qlen) {
#             Get scores in FASTQ format (Phred+33)
              Q = qualS[++qpos]
              qline = qline sprintf("%c", (Q <= HIGH_CHAR ? Q : HIGH_CHAR) + LOW_CHAR)
            }
            qualA[readID++] = sprintf("+\n%s\n", qline)
            qual = ""
            qline = ""
          }
        }
      }
      if (($0 ~ /^}/) && (runId != prevRunId)) {
        for (pairNum = 1; pairNum <= readID; pairNum++) {
          printf("%s", pairA[pairNum]) >> outFile
          printf("%s", qualA[pairNum]) >> outFile
        }
#       This prevent duplication last Read
        prevRunId = runId
      }
    }
    close(query)
    print "[" p "/" l "] recd " pageSize " spots from " firstRID " to " runId " on " currPage " page"
    system("sleep 8")
  }
  close(outFile);
  print "\nWritten " j " sequences from " pages " pages at " runAcc
}

function ceil(n,    _in)
{
  _in = int(n)
  if (n == _in) {
    return _in
  }
  return _in + 1
}

function min(x, y)
{
  return (x < y) ? x : y
}

function gcd(a, b)
{
    while (a != b)
        if (a > b)
           a = a - b
        else
           b = b - a
    return a
}
