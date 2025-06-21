<script setup>
/* eslint-disable */
import { ref, onMounted, watch } from "vue";
import { useRoute } from "vue-router";
import { backendLink } from "@/scripts/AccountFunc";
import axios from "axios";
import histogramData from "../data/histogram.json";
import sunburstTableData from "../data/sunburst.table.json";
import * as echarts from "echarts";
import "echarts/lib/chart/sunburst";

const route = useRoute();
const sunburstRef = ref(null);

//// HISTOGRAM ////
const lociDropdownValue = ref("HLA-B");
const drugDropdownValue = ref("carbamazepine");

const lociDropdown = ref();
const drugDropdown = ref();
const regionDropdown = ref(); // List of region names

const selectedRegions = ref([false]); // Checkbox v-model
const selectedRegionsRef = ref([]); // List of ref of charts


const histogramRefMany = ref();
const histogramRefThai = ref();
//// HISTOGRAM ////

//// SUNBURST ////
let sunburstData;
const minPRef = ref(false);
const sunburstSelectedCountries = ref([]); // Checkbox v-model
const sunburstCountryDropdown = ref([]); // List of country names

const sunburstFilter = ref([]);
const sunburstTableContent = ref([]);
//// SUNBURST ////

async function GetSunburstData() {
  /*
  Object(
    country_filter:Array<String>,
    minP_filter:Bool
    )
  */
  let selectedCountries = sunburstCountryDropdown.value.filter(
    (item, index) => sunburstSelectedCountries.value[index + 1]
  );
  axios.post(backendLink + "/data/sunburst", {
    country_filter: selectedCountries,
    minP_filter: minPRef.value,
  }).then((response) => {
    sunburstData = response.data;
    SunburstChart();
    SunburstTable();
  }).catch((error) => console.error(error));
  
}

function HistogramFilterSetup() {
  const lociSet = new Set();
  const drugSet = new Set();
  const regionSet = new Set();
  for (let obj of histogramData) {
    lociSet.add(obj.loci);
    for (let drug of obj.drug) if (drug != "") drugSet.add(drug);

    regionSet.add(obj.region);
  }

  lociDropdown.value = Array.from(lociSet);
  drugDropdown.value = Array.from(drugSet);
  regionDropdown.value = Array.from(regionSet);
  for (let i = 0; i < regionDropdown.value.length; i++) {
    selectedRegions.value.push(false);
  }
}

function SunburstFilterSetup() {
  const regionSet = new Set();
  for (let obj of sunburstTableData) regionSet.add(obj.country);
  sunburstCountryDropdown.value = Array.from(regionSet);
  sunburstSelectedCountries.value = Array(
    sunburstCountryDropdown.value.length + 1
  ).fill(true);
}

function SunburstTable(params) {
  if (params && params.data) {
    let obj = params.data;
    if (!obj.name) sunburstFilter.value.pop();
    else {
      sunburstFilter.value = structuredClone(obj.parent);
      sunburstFilter.value.push(obj.name);
    }
  }

  sunburstTableContent.value = [];
  for (let obj of sunburstTableData) {
    // Filter for the sunburst chart click
    let bool_value = true;
    const levels = ["disease", "drug", "adr", "gene_loci", "allele"];
    for (let index in sunburstFilter.value) {
      bool_value =
        bool_value && sunburstFilter.value[index] == obj[levels[index]];
    }
    if (!bool_value) continue;

    // Filter for the country checkbox
    bool_value = false;
    for (let index in sunburstSelectedCountries.value) {
      if (index == 0) continue;
      if (
        sunburstSelectedCountries.value[index] &&
        sunburstCountryDropdown.value[index - 1] == obj.country
      ) {
        bool_value = true;
        break;
      }
    }
    if (!bool_value) continue;

    // Filter for MinP
    if (minPRef.value && !obj.minP) continue;

    sunburstTableContent.value.push([
      obj.drug,
      obj.allele,
      obj.adr,
      obj.ethnicity,
      obj.pvalue_exposed,
      obj.pvalue_population,
      obj.odds_exposed,
      obj.pubmed,
      obj.allele_count,
      obj.allele_freq,
    ]);
  }
}

async function SunburstChart() {
  if (echarts.getInstanceByDom(sunburstRef.value)) echarts.dispose(echarts.getInstanceByDom(sunburstRef.value));
  const chart = echarts.init(sunburstRef.value);
  const option = {
    tooltip: {
      trigger: "item",
      formatter: "{b}: {c}",
    },
    series: {
      type: "sunburst",
      itemStyle: { color: "#242756" },
      data: sunburstData,
      radius: [0, "90%"],
      label: {
        show: false,
      },
      levels: [
        { itemStyle: { color: "#242756" } },
        { itemStyle: { color: "#DB6F57" } },
        { itemStyle: { color: "#ECE0C5" } },
        { itemStyle: { color: "#CF648C" } },
        { itemStyle: { color: "#5DA597" } },
        { itemStyle: { color: "#55547F" } },
      ],
    },
  };
  chart.setOption(option);

  chart.on("click", SunburstTable);
  window.addEventListener("resize", chart.resize);
}

function HistogramChart2() {
  const initChart = (histogramElement, region) => {
    if (echarts.getInstanceByDom(histogramElement))
      echarts.dispose(echarts.getInstanceByDom(histogramElement));

    const filteredData = histogramData.filter((item) => {
      let asdf = true;
      if (drugDropdownValue.value != "All")
        asdf = asdf && item.drug.includes(drugDropdownValue.value);
      if (region != "All") asdf = asdf && item.region == region;
      if (lociDropdownValue.value != "All")
        asdf = asdf && item.loci == lociDropdownValue.value;
      return asdf;
    });
    const chart = echarts.init(histogramElement);
    if (filteredData.length == 0) {
      chart.clear();
      return;
    }
    const sortedData = filteredData.sort((a, b) => b.data[0] - a.data[0]);

    const option = {
      tooltip: {
        trigger: "axis",
        axisPointer: {
          type: "shadow-sm",
        },
      },
      legend: {},
      xAxis: {
        type: "value",
        min: 0, // Set the minimum value for the scale
        max: 100, // Set the maximum value for the scale
        interval: 10, // Define intervals for the scale
      },
      yAxis: {
        type: "category",
        data: [`${region} (%)`], // Single bar
      },
      series: sortedData,
    };

    chart.setOption(option);
  };


  let thaiIndex = regionDropdown.value.findIndex((element) => element == "TH_AF")
  if (selectedRegions.value[thaiIndex+1]) {
    initChart(
      histogramRefThai.value,
      "TH_AF"
    );
  }
  
}

function HistogramChart() {
  var result = [];
  var all_type = [];
  var count_stack = [];
  if (echarts.getInstanceByDom(histogramRefMany.value))
  echarts.dispose(echarts.getInstanceByDom(histogramRefMany.value));
  const initChart = (region) => {
    if (region != "TH_AF") {
      const filteredData = histogramData.filter((item) => {
        let asdf = true;
        if (drugDropdownValue.value != "All")
          asdf = asdf && item.drug.includes(drugDropdownValue.value);
        if (region != "All") asdf = asdf && item.region == region;
        if (lociDropdownValue.value != "All")
          asdf = asdf && item.loci == lociDropdownValue.value;
        return asdf;
      });
      if (filteredData.length == 0) return;
      
      const sortedData = filteredData.sort((a, b) => b.data[0] - a.data[0]);

      var reg = { product: sortedData[0]["region"] };
      sortedData.forEach((element) => {
        reg[element.name] = element.data[0];
        all_type.push(element.name);
      });
      result.push(reg);
    }
  };
  for (let asdf in selectedRegions.value) {
    if (selectedRegions.value[asdf] && asdf != 0) {
      initChart(
        regionDropdown.value[asdf - 1]
      );
    }
  }
  // console.log(result);
  all_type = [...new Set(all_type)];
  // console.log(all_type);
  all_type.forEach((element) => {
    count_stack.push({ type: "bar", stack: "q" });
  });
  all_type.unshift("product");

  var myChart = echarts.init(histogramRefMany.value);
  var option;

  option = {
    legend: {},
    tooltip: {
      trigger: "axis",
      axisPointer: {
        type: "shadow-sm",
      },
    },
    dataset: {
      dimensions: all_type,
      source: result,
    },
    xAxis: {
      type: "value",
      min: 0, // Set the minimum value for the scale
      max: 100, // Set the maximum value for the scale
      interval: 10, // Define intervals for the scale
    },
    yAxis: { type: "category" ,
    },
    series: count_stack,
  };
  myChart.setOption(option);
}

function SelectRegion(region) {
  if (region == "All") {
    let value = !selectedRegions.value[0];
    for (let i in selectedRegions.value) {
      selectedRegions.value[i] = value;
    }
  } else if (region == "TH_AF") {
    let index = regionDropdown.value.findIndex((asdf) => asdf == "TH_AF");
    selectedRegions.value[index + 1] = true;
  }

  HistogramChart();
  HistogramChart2();
}

function SelectCountry(country) {
  if (country == "All") {
    let value = !sunburstSelectedCountries.value[0];
    for (let i in sunburstSelectedCountries.value) {
      sunburstSelectedCountries.value[i] = value;
    }
  } else if (country == "Thailand") {
    let index = sunburstCountryDropdown.value.findIndex(
      (asdf) => asdf == "Thailand"
    );
    sunburstSelectedCountries.value[index + 1] = true;
  }
  GetSunburstData();
}

function InitRef(el) {
  if (selectedRegionsRef.value.length < regionDropdown.value.length)
    selectedRegionsRef.value.push(el);
}

function ActivateSearch() {
  if (route.query.q) {
    drugDropdownValue.value = route.query.q;
    lociDropdownValue.value = "All";
  }
}

onMounted(ActivateSearch);

watch(() => route.fullPath, ActivateSearch);
onMounted(HistogramFilterSetup);
onMounted(SunburstFilterSetup);
onMounted(SunburstTable);
onMounted(SunburstChart);
onMounted(() => setTimeout(() => SelectRegion("All"), 300));
onMounted(GetSunburstData);

function WELOVEVUEJSBABY(link) {
  window.location.href = link;
}
</script>

<template>
  <div class="absolute top-20 grid gap-16 grid-cols-1 w-full py-10">
    <div class="flex justify-center gap-10 items-center w-full">
      <div
        class="w-11/12 max-w-300 min-w-96 h-full shadow-md bg-white rounded-md transition-all duration-500 ease-in-out"
      >
        <div
          class="h-20 w-full justify-center items-center flex font-bold text-2xl"
        >
          Thailand Allele Frequency Chart
        </div>
        <div class="w-full h-64" ref='histogramRefThai'></div>
        <div class="w-full h-144" ref='histogramRefMany'></div>
        <div class="w-full h-16 flex justify-center items-center gap-10">
          <details class="w-64 h-14 flex justify-start dropdown">
            <summary class="btn bg-white border-primary">
              Select regions
            </summary>
            <ul
              multiple
              class="menu dropdown-content bg-white rounded-box z-1 p-2 shadow-sm w-200 h-32"
            >
              <li>
                <label class="cursor-pointer justify-between flex items-center">
                  <span> All </span>
                  <input
                    type="checkbox"
                    checked="checked"
                    @change="SelectRegion('All')"
                    class="checkbox checkbox-primary"
                  />
                </label>
              </li>
              <li v-for="(i, index) in regionDropdown" :key="i">
                <label class="cursor-pointer flex justify-between items-center">
                  <span class="">{{ i }}</span>
                  <input
                    type="checkbox"
                    v-model="selectedRegions[index + 1]"
                    @change="SelectRegion()"
                    class="checkbox checkbox-primary"
                  />
                </label>
              </li>
            </ul>
          </details>
          <div class="w-64 h-14 flex justify-start">
            <select
              v-model="drugDropdownValue"
              @change="SelectRegion"
              class="select select-primary w-full max-w-xs"
            >
              <option value="All">All</option>
              <option v-for="i in drugDropdown" :key="i" :value="i">
                {{ i }}
              </option>
            </select>
          </div>
          <div class="w-64 h-14 flex justify-start">
            <select
              v-model="lociDropdownValue"
              @change="SelectRegion"
              class="select select-primary w-full max-w-xs"
            >
              <option value="All">All</option>
              <option v-for="i in lociDropdown" :key="i" :value="i">
                {{ i }}
              </option>
            </select>
          </div>
        </div>
      </div>
    </div>

    <div class="flex justify-center items-center w-full">
      <div class="w-11/12 max-w-360 max-h-512 shadow-md bg-white rounded-md">
        <h1
          class="w-full h-32 text-4xl font-bold text-center flex justify-center items-center"
        >
          TH~DrugXcape
        </h1>

        <div class="w-full h-220 flex justify-center">
          <div ref="sunburstRef" class="w-220 h-220"></div>
          <div class="w-64 h-full grid grid-cols-1">
            <div class="w-full h-160 overflow-y-scroll overflow-x-hidden">
              <details open class="w-64 flex justify-start dropdown dropdown-bottom">
                <summary class="hidden">
                  Select country
                </summary>
                <ul
                  multiple
                  class="w-full menu dropdown-content dropdown-bottom bg-white overflow-y-scroll rounded-box z-1 p-2 shadow-sm"
                >
                  <li>
                    <label class="flex justify-between items-center cursor-pointer">
                      <span class="label-text"> All </span>
                      <input
                        type="checkbox"
                        checked="checked"
                        @change="SelectCountry('All')"
                        class="checkbox checkbox-primary"
                      />
                    </label>
                  </li>
                  <li v-for="(i, index) in sunburstCountryDropdown" :key="i">
                    <label class="flex justify-between items-center cursor-pointer">
                      <span class="label-text">{{ i }}</span>
                      <input
                        type="checkbox"
                        v-model="sunburstSelectedCountries[index + 1]"
                        @change="SelectCountry()"
                        class="checkbox checkbox-primary"
                      />
                    </label>
                  </li>
                </ul>
              </details>
            </div>
            <label
              class="flex justify-around items-center cursor-pointer w-64 h-16 border-2 border-accent"
            >
              <span class="label-text pl-2">MinP</span>
              <input
                type="checkbox"
                class="checkbox checkbox-accent"
                v-model="minPRef"
                @change="SelectCountry()"
              />
            </label>
          </div>
        </div>
        <div class="w-full h-8 flex justify-center items-center">
          <div class="w-220 h-full justify-center flex items-center gap-10">
            <div class="flex justify-center items-center gap-1">
              <div class="w-5 h-3 bg-[#DB6F57] rounded-xs"></div>
              <div class="text-sm font-semibold">Disease</div>
            </div>
            <div class="flex justify-center items-center gap-1">
              <div class="w-5 h-3 bg-[#ECE0C5] rounded-xs"></div>
              <div class="text-sm font-semibold">Drug</div>
            </div>
            <div class="flex justify-center items-center gap-1">
              <div class="w-5 h-3 bg-[#CF648C] rounded-xs"></div>
              <div class="text-sm font-semibold">ADR</div>
            </div>
            <div class="flex justify-center items-center gap-1">
              <div class="w-5 h-3 bg-[#5DA597] rounded-xs"></div>
              <div class="text-sm font-semibold">Gene Loci</div>
            </div>
            <div class="flex justify-center items-center gap-1">
              <div class="w-5 h-3 bg-[#55547F] rounded-xs"></div>
              <div class="text-sm font-semibold">Allele</div>
            </div>
          </div>
          <div class="w-64"></div>
        </div>

        <div class="py-2"></div>
        <div class="w-full max-h-200 overflow-x-auto">
          <table
            class="table table-sm table-pin-rows table-pin-cols text-center overflow-scroll h-full"
          >
            <thead>
              <tr>
                <th></th>
                <td>Drug</td>
                <td>Allele</td>
                <td>ADR</td>
                <td>Ethnicity</td>
                <td>PValue Exposed</td>
                <td>PValue Population</td>
                <td>Odds Exposed</td>
                <td>Publication (Pubmed)</td>
                <td>Allele Count</td>
                <td>Allele Frequency (%)</td>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(row, index) in sunburstTableContent" :key="index">
                <th>{{ index + 1 }}</th>
                <td>{{ row[0] }}</td>
                <td>{{ row[1] }}</td>
                <td>{{ row[2] }}</td>
                <td>{{ row[3] }}</td>
                <td>{{ row[4] }}</td>
                <td>{{ row[5] }}</td>
                <td>{{ row[6] }}</td>
                <td @click="WELOVEVUEJSBABY(row[7])" class="link">
                  {{ row[7] }}
                </td>
                <td>{{ row[8] }}</td>
                <td>{{ row[9] }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- Search Results 
  <div v-if="route.query.q" class="flex justify-center items-center min-h-96 w-full">
    <div class="w-11/12 min-w-96 max-w-256 h-full shadow-md bg-white rounded-md">
      <div class="w-full h-20 flex items-center justify-normal">
        <h1 class="px-8 pt-4 text-3xl">You searched for : <span class="font-bold">{{ route.query.q }}</span></h1>
      </div>
      <div class="divider h-0"></div>
      <div v-for="i in [1,2,3]" :key="i" class=" h-64 py-2 flex justify-center items-center">
        <div class="w-49/50 h-full shadow-xl rounded-xl bg-white">
          <h1 class="text-lg w-full px-6 py-3">DRUG {{ i }}</h1>
        </div>
      </div>
    </div>  
  </div>
  Search Results -->
  </div>
</template>
