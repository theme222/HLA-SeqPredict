/* eslint-disable */


function errorHandler(error)
{
    window.sharedData.latestError = error
    console.error(error)
}

let igvMade = 0
let currentBrowser
let currentTrack 
let igvDivDiv = document.getElementById("igv-div-div")
let igvDiv = document.getElementById("igv-div");
let prevUrl = null
let prevSequence = ""
let backendLink = window.sharedData.backendLink // fuckin lifesaver

async function makeIGV(igvDiv) { 
    if (igvMade == 1) return;
    igvMade += 1
    var options = {
        reference: {
            "id": "chr6",
            "name": "Human Chromosome 6",
            "fastaURL": `${backendLink}/reference/chr6.fa`,
            'indexURL': `${backendLink}/reference/chr6.fa.fai`,
        },
        locus: 'NC_000006.12:29,940,368-29,947,527'
    };

    igv.createBrowser(igvDiv, options).then(function (browser) {
        currentBrowser = browser
        console.log("Created IGV browser");
        console.log(igvMade)
    }).catch(errorHandler);
}

function interval()
{
    // Displaying igv
    let currURL = window.location.href
    let onDashboard = window.location.href.split('/').pop() == 'dashboard'
    if (currURL != prevUrl) {
        prevUrl = currURL
        // Checking if on dashboard screen then hide if not
        if (onDashboard) { 
            igvDivDiv.style.display = 'block'
        } else {
            igvDivDiv.style.display = 'none'
        }
    }

    // Getting current sequence
    if (!onDashboard) {return}
    let currentSequence = document.getElementById("selectedSequence").textContent
    if (currentSequence == "") return
    if (currentSequence != prevSequence && currentBrowser)
    {
        prevSequence = currentSequence
        axios.post(`${backendLink}/api/requestFile/igv`, // TODO: MAKE SURE TO ADD /igv AT THE END AS WELLL
        {
            session_cookie: Cookies.get('session_cookie'),
            label: currentSequence
        })
        .then(response => {
            console.log(`loading track info from ${response.data['token_bam']} and ${response.data['token_bam_bai']}`)
            
            if (currentTrack){
                currentBrowser.removeTrack(currentTrack)
            }
            console.log(response.data['range'])
            currentBrowser.search("NC_000006.12:"+response.data['range'])
            currentBrowser.loadTrack(
                {
                    "name": "Your sequence",
                    "format": "bam",
                    "url": `${backendLink}/download/` + response.data['token_bam'],
                    "indexURL": `${backendLink}/download/` + response.data['token_bam_bai']
                    })
                .then(function (newTrack) {
                alert("Track load success")
                currentTrack = newTrack
                })
        }).catch(errorHandler)
    }


}

makeIGV(igvDiv)
setInterval(interval, 100)

