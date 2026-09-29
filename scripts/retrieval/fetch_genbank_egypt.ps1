<#
  fetch_ncbi_ALL.ps1
  Comprehensive NCBI GenBank retrieval for Egyptian H5 and H9 avian influenza.

  What is different from the earlier fetch_ncbi_HA.ps1:
    * NO length filter. Every record is retrieved, partial and complete.
      Filtering happens downstream where it can be documented.
    * STRUCTURED METADATA for every accession - country, host, collection date,
      serotype, segment, length - written to a TSV beside the FASTA. This is the
      field that was missing last time. Country from the record's own source
      qualifier is checkable; country guessed from a strain name is not.
    * NA (segment 6) as well as HA (segment 4), so the N8 -> N1 -> N2 shift can be
      confirmed on the neuraminidase sequence itself rather than read off the
      subtype string in the definition line.

  Run it from the folder you want the files in:
      cd "C:\Users\ahme_\Claude\Projects\LLM"
      .\fetch_ncbi_ALL.ps1
  If PowerShell blocks it:
      powershell -ExecutionPolicy Bypass -File .\fetch_ncbi_ALL.ps1

  Outputs:
      NCBI_H5_HA.fasta   NCBI_H5_HA_metadata.tsv
      NCBI_H9_HA.fasta   NCBI_H9_HA_metadata.tsv
      NCBI_H5_NA.fasta   NCBI_H5_NA_metadata.tsv     (if $FetchNA)
      NCBI_H9_NA.fasta   NCBI_H9_NA_metadata.tsv     (if $FetchNA)
      NCBI_fetch_log.txt
#>

# =========================================================================
#  RUN THIS AS A FILE, NOT BY PASTING IT INTO THE CONSOLE.
#      cd "C:\Users\ahme_\Claude\Projects\LLM"
#      .\fetch_ncbi_ALL.ps1
#  Pasting the body into the prompt runs it in whatever folder you happen to
#  be in and mangles the multi-line blocks.
# =========================================================================

# ----------------------------- SETTINGS -----------------------------------
$Email   = "your.email@example.org"   # NCBI asks for a contact address on bulk queries.
$ApiKey  = ""        # <-- optional NCBI API key; raises the rate limit 3/sec -> 10/sec.
                     #     Free from https://account.ncbi.nlm.nih.gov/settings/
$FetchNA = $true     # $false to skip the neuraminidase segments
$Batch   = 200       # records per request
# --------------------------------------------------------------------------

$ErrorActionPreference = "Stop"

# Windows PowerShell 5.1 still defaults to TLS 1.0 on some systems; NCBI refuses it.
[Net.ServicePointManager]::SecurityProtocol =
    [Net.SecurityProtocolType]::Tls12 -bor [Net.SecurityProtocolType]::Tls11

$base    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Write output beside the script itself, so the folder you launched from
# does not matter. Falls back to the current folder if pasted interactively.
if ($PSScriptRoot) { $here = $PSScriptRoot } else { $here = (Get-Location).Path }
Write-Host "Output folder: $here" -ForegroundColor Cyan
$log     = Join-Path $here "NCBI_fetch_log.txt"
$utf8    = New-Object System.Text.UTF8Encoding($false)   # no BOM - BOM breaks TSV readers
"=== NCBI fetch started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
    Out-File $log -Encoding ascii

function Write-Log($msg) {
    Write-Host $msg
    $msg | Out-File $log -Append -Encoding ascii
}

$extra = "&tool=EgyptAIV_survey"
if ($Email)  { $extra += "&email=$([uri]::EscapeDataString($Email))" }
if ($ApiKey) { $extra += "&api_key=$ApiKey" }
$delayMs = if ($ApiKey) { 110 } else { 350 }

# Read one direct child element's text, or "" if absent.
function Get-NodeText($node, $name) {
    if (-not $node) { return "" }
    $n = $node.SelectSingleNode($name)
    if ($n) { return $n.InnerText } else { return "" }
}

function Invoke-Eutil($url) {
    for ($try = 1; $try -le 4; $try++) {
        try {
            Start-Sleep -Milliseconds $delayMs
            return Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 300
        } catch {
            if ($try -eq 4) { throw }
            $wait = 5 * $try
            Write-Log "    request failed (attempt $try of 4): $($_.Exception.Message)"
            Write-Log "    retrying in $wait s ..."
            Start-Sleep -Seconds $wait
        }
    }
}

function Build-Query($subtypes, $segmentNo, $geneNames) {
    $sub = ($subtypes  | ForEach-Object { "`"$_`"[All Fields]" }) -join " OR "
    $gen = ($geneNames | ForEach-Object { "`"$_`"[All Fields]" }) -join " OR "
    return "(`"Influenza A virus`"[Organism]) AND (Egypt[All Fields]) AND ($sub) AND (($gen) OR (`"segment $segmentNo`"[All Fields]))"
}

function Fetch-Set($label, $query, $fastaOut, $metaOut) {
    Write-Log ""
    Write-Log "=== $label ==="
    Write-Log "  query: $query"

    # 1. esearch onto the history server
    $u = "$base/esearch.fcgi?db=nuccore&term=$([uri]::EscapeDataString($query))&usehistory=y&retmax=0$extra"
    # PowerShell's XML adapter can flatten $xml.eSearchResult to a String, so use
    # the underlying .NET XmlDocument API rather than dotted property access.
    $doc = New-Object System.Xml.XmlDocument
    $doc.LoadXml((Invoke-Eutil $u).Content)
    $root   = $doc.DocumentElement                 # <eSearchResult>
    $count  = [int](Get-NodeText $root "Count")    # direct child, not the one in TranslationStack
    $webenv = Get-NodeText $root "WebEnv"
    $qkey   = Get-NodeText $root "QueryKey"
    if (-not $webenv -or -not $qkey) { throw "esearch did not return a history handle; check the query." }
    Write-Log "  $count records found"

    $header = "accession`tlength`tcountry`thost`tcollection_date`tserotype`tsegment`tstrain`ttitle"
    if ($count -eq 0) {
        [System.IO.File]::WriteAllText($fastaOut, "", $utf8)
        [System.IO.File]::WriteAllLines($metaOut, @($header), $utf8)
        return 0
    }

    # 2. efetch FASTA in batches
    $sw = New-Object System.IO.StreamWriter($fastaOut, $false, [System.Text.Encoding]::ASCII)
    try {
        for ($i = 0; $i -lt $count; $i += $Batch) {
            $u = "$base/efetch.fcgi?db=nuccore&query_key=$qkey&WebEnv=$webenv&retstart=$i&retmax=$Batch&rettype=fasta&retmode=text$extra"
            $sw.Write((Invoke-Eutil $u).Content)
            Write-Log ("    sequences {0}/{1}" -f [Math]::Min($i + $Batch, $count), $count)
        }
    } finally { $sw.Close() }

    # 3. esummary in batches -> structured source qualifiers.
    #    SubType holds the qualifier names, SubName the matching values, pipe-separated.
    $rows = New-Object System.Collections.Generic.List[string]
    $rows.Add($header)
    for ($i = 0; $i -lt $count; $i += $Batch) {
        $u = "$base/esummary.fcgi?db=nuccore&query_key=$qkey&WebEnv=$webenv&retstart=$i&retmax=$Batch&version=2.0$extra"
        $sdoc = New-Object System.Xml.XmlDocument
        $sdoc.LoadXml((Invoke-Eutil $u).Content)
        foreach ($d in $sdoc.GetElementsByTagName("DocumentSummary")) {
            if (-not $d) { continue }
            $acc   = Get-NodeText $d "AccessionVersion"
            if (-not $acc) { $acc = Get-NodeText $d "Caption" }
            $slen  = Get-NodeText $d "Slen"
            $title = (Get-NodeText $d "Title") -replace "[`t`r`n]", " "
            $q = @{}
            $st = Get-NodeText $d "SubType"
            $sn = Get-NodeText $d "SubName"
            if ($st -and $sn) {
                $keys = $st -split "\|"
                $vals = $sn -split "\|"
                for ($k = 0; $k -lt [Math]::Min($keys.Count, $vals.Count); $k++) {
                    if ($keys[$k]) { $q[$keys[$k]] = ($vals[$k] -replace "[`t`r`n]", " ") }
                }
            }
            $strain = ""
            if     ($title -match "\(([^()]*?)\s*\(H\d+N\d+\)\)") { $strain = $Matches[1] }
            elseif ($title -match "\((A[/\-][^()]+)\)")           { $strain = $Matches[1] }
            $get = { param($k) if ($q.ContainsKey($k)) { $q[$k] } else { "" } }
            $rows.Add( ($acc, $slen,
                        (& $get "country"), (& $get "host"),
                        (& $get "collection_date"), (& $get "serotype"),
                        (& $get "segment"), $strain, $title) -join "`t" )
        }
        Write-Log ("    metadata  {0}/{1}" -f [Math]::Min($i + $Batch, $count), $count)
    }
    [System.IO.File]::WriteAllLines($metaOut, $rows, $utf8)

    $n = @(Select-String -Path $fastaOut -Pattern "^>").Count
    Write-Log "  wrote $n sequences    -> $(Split-Path $fastaOut -Leaf)"
    Write-Log "  wrote $($rows.Count - 1) metadata rows -> $(Split-Path $metaOut -Leaf)"
    if ($n -ne ($rows.Count - 1)) {
        Write-Log "  NOTE: sequence count and metadata count differ. Join on accession downstream."
    }
    return $n
}

$H5 = @("H5N1","H5N2","H5N3","H5N4","H5N5","H5N6","H5N7","H5N8","H5N9","H5Nx","H5")
$H9 = @("H9N2","H9N1","H9N3","H9N6","H9N9","H9Nx","H9")
$HA = @("hemagglutinin","haemagglutinin","HA")
$NA = @("neuraminidase","NA")

$total = 0
$total += Fetch-Set "H5 HA (segment 4), Egypt, no length filter" (Build-Query $H5 4 $HA) `
            (Join-Path $here "NCBI_H5_HA.fasta") (Join-Path $here "NCBI_H5_HA_metadata.tsv")
$total += Fetch-Set "H9 HA (segment 4), Egypt, no length filter" (Build-Query $H9 4 $HA) `
            (Join-Path $here "NCBI_H9_HA.fasta") (Join-Path $here "NCBI_H9_HA_metadata.tsv")
if ($FetchNA) {
    $total += Fetch-Set "H5 NA (segment 6), Egypt" (Build-Query $H5 6 $NA) `
                (Join-Path $here "NCBI_H5_NA.fasta") (Join-Path $here "NCBI_H5_NA_metadata.tsv")
    $total += Fetch-Set "H9 NA (segment 6), Egypt" (Build-Query $H9 6 $NA) `
                (Join-Path $here "NCBI_H9_NA.fasta") (Join-Path $here "NCBI_H9_NA_metadata.tsv")
}

Write-Log ""
Write-Log "=== done: $total sequences across all sets ==="
Write-Log "Files are in $here"
Write-Host ""
Write-Host "Finished. Tell Claude the files are ready." -ForegroundColor Green
