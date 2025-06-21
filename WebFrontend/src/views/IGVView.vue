<script setup>
/* eslint-disable */
import { backendLink } from '@/scripts/AccountFunc';
import { onMounted, ref } from 'vue';
import { RouterLink, useRoute } from 'vue-router';
import Cookies from 'js-cookie';
import axios from "axios";
import { ErrorHandler } from '@/scripts/OtherFunc';

const route = useRoute();

const sequence_id = route.query.sequence_id;
const igvDiv = ref(null);
const showAlignments = ref(false);

let currentBrowser;
let currentTrack;


const browserOptions = {
    reference: {
        "id": "chr6",
        "name": "Human Chromosome 6",
        "fastaURL": `${backendLink}/reference/chr6.fa`,
        'indexURL': `${backendLink}/reference/chr6.fa.fai`,
    },
    locus: 'NO_000006.12:29,940,368-29,947,527'
};

async function MakeIGVBrowser()
{
    currentBrowser = await igv.createBrowser(igvDiv.value, browserOptions).catch(ErrorHandler);
    if (!currentBrowser) return;

    console.log("Created IGV browser");
    LoadSequence();
}

async function LoadSequence()
{
    if (!currentBrowser) return;
    if (!sequence_id) return;

    let response = await axios.post(`${backendLink}/api/requestFile/igv`,
        {
            cookie: Cookies.get("session_cookie"),
            sequence_id: sequence_id,
        }
    ).catch(ErrorHandler);

    if (!response) return;
    let data = response.data;

    if (currentTrack) currentBrowser.removeTrack(currentTrack);

    let track = {
        "name": "Your sequence",
        "format": "bam",
        "url": `${backendLink}/download/` + data['token_bam'],
        "indexURL": `${backendLink}/download/` + data['token_bam_bai'],
        "showAlignments": showAlignments.value, 
    }

    currentTrack = await currentBrowser.loadTrack(track).catch(ErrorHandler);
    currentBrowser.search("NC_000006.12:" + data['range']);
    if (!currentTrack) return;
}

onMounted(MakeIGVBrowser);

</script>
<template>

  <div class="w-full absolute top-30 grid grid-cols-1 gap-10">
    <div class="w-full flex justify-center items-center">
        <div class="w-11/12 bg-white rounded-md h-150 shadow-md" id="infoPage">
            <div class="w-full h-30 flex gap-10 justify-center items-center">
                <h1 class="text-4xl font-bold py-6 h-30 flex items-center justify-center">Visualize files with&nbsp;<a href="https://igv.org/" class="link">IGV.js</a></h1>
                <button @click="LoadSequence" class="btn btn-dash btn-primary">Reload</button>
            </div>
            <div class="w-full h-120" ref="igvDiv" id="igvDiv"></div>
        </div>
    </div>  

    <div class="w-full h-20 flex justify-center items-center">
        <div class="w-11/12 h-full bg-white rounded-md shadow-xl flex justify-center items-center">
            <div class="w-48/50 flex gap-6 items-center h-full">
                <h1 class="text-2xl font-semibold">Options :</h1> 
                <label class="w-50 h-16 rounded-md shadow-sm flex justify-center items-center gap-4 bg-base-300 cursor-pointer">
                    <span>Show Alignments</span>
                    <input type="checkbox" class="checkbox checkbox-primary" id="showAlignmentCheckbox" v-model="showAlignments">
                </label>
            </div>
        </div>
    </div>

  </div>


</template>

