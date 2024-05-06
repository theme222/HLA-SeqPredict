/* eslint-disable */
let igvMade = 0
let currentBrowser
let currentTrack 


async function makeIGV(igvDiv) { /* 
    let seqURL = await axios.post("http://localhost:7000/api/getfile")
                .catch(error => {
                    alert(error)
                    console.log("Error getting file")
                })

    console.log(seqURL)
    */
    var options = {
        reference: {
            "id": "chr6",
            "name": "Human Chromosome 6",
            "fastaURL": "http://127.0.0.1:7000/reference/chr6.fa",
            'indexURL': "http://127.0.0.1:7000/reference/chr6.fa.fai",

            /* 
            "tracks": [
                {
                    "name": "Your sequence",
                    "format": "bam",
                    "url": "http://127.0.0.1:7000/api/getfile",
                    "indexURL": "http://127.0.0.1:7000/api/getfileaaa"
                }
            ]
            */
        },
        locus: 'NC_000006.12:29,940,368-29,947,527'
    };

    igv.createBrowser(igvDiv, options).then(function (browser) {
        currentBrowser = browser
        console.log("Created IGV browser");
        igvMade += 1
    }).catch(function (error) {
        console.error("Error creating IGV browser:", error);
    });
}

let prevUrl = null
setInterval(() => {
    let currURL = window.location.href
    if (currURL != prevUrl) {
        prevUrl = currURL
        if (window.location.href.split('/').pop() == 'dashboard') { // last element in array is dashboard it is the correct page

            var igvDiv = document.getElementById("igv-div");
            igvDiv.style.zIndex = 1
            if (igvMade == 0) {
                makeIGV(igvDiv)
            }
            document.getElementById("visualizeIGV").onclick = function () {
                if (currentBrowser) {
                    let currentSequence = document.getElementById("sequences").value
                    if (currentSequence) {
                        axios.post("http://localhost:7000/api/requestFile", {
                            session_cookie: Cookies.get('session_cookie'),
                            label: currentSequence
                        }).then(response => {
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
                                    "url": "http://localhost:7000/download/" + response.data['token_bam'],
                                    "indexURL": "http://localhost:7000/download/" + response.data['token_bam_bai']
                                  })
                              .then(function (newTrack) {
                                alert("Track load success")
                                currentTrack = newTrack
                              })
                              .catch(function (error)  {
                                 alert(error)
                                 console.error(error)
                              })

                        }).catch(error => {
                            alert(error.data)
                            console.log(error.data)
                        })
                    }

                }
            };


        } else {
            let igvDiv = document.getElementById("igv-div");
            igvDiv.style.zIndex = -1
        }
    }
}, 1)
