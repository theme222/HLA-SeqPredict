<template>
  <div></div>
</template>

<script setup>
/* eslint-disable */
import pdfMake from "pdfmake"; // main pdf lib
import router from "@/router";
import { sharedData } from "@/scripts/SharedData";
import { text } from "@fortawesome/fontawesome-svg-core";
// import pdfFonts from "@/assets/font/custom-fonts.js"; // add font

/*
[
      {
        text: [
          { text: "Drug: ", style: "titleBold" },
          { text: "[Drug Name]" },
        ],
      },
      {
        ul: [
          {
            text: "ADR: [Specific ADR, e.g., Stevens-Johnson Syndrome]",
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          },
          {
            text: "HLA Allele Association: [Specific HLA allele(s) associated with the ADR]",
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          },
          {
            text: "Quality Metric: [Value of quality metric]",
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          },
        ],
      },
    ],
*/

function formatADR()
{
  let returnList = []
  for (let item of sharedData.drugList)
  {
    returnList.push(
      {
        text: [
          { text: "Drug: ", style: "titleBold" },
          { text: `${item.drug} : ${item.pubmed}` },
        ],
      }
    )

    let textObject = {
        ul: [
          {
            text: `ADR: ${item.adr}`,
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          },
          {
            text: `Associated Allele: ${item.allele}`,
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          },
        ],
      }

      if (item.odds_exposed)
      {
        textObject.ul.push({
            text: `Odds: ${item.odds_exposed}`,
            listType: "circle",
            margin: [20, 0, 0, 0],
            style: "normalText",
          })
      }    

    returnList.push(textObject)
  }
  return returnList
}

async function pdfRender() {
      let docDefinition = {
        pageSize: "A4",
        // pageMargins: [30, 20, 30, 90],
        content: [
          {
            text: "HLA Genotype-Based Adverse Drug Reaction Prediction",
            style: "header",
          },
          {text:" "},
          {
            text: [
              { text: "Sequence Label: ", style: "titleBold" },
              { text: sharedData.sequenceLabel, fontSize: 11 },
            ],
          },
          {
            text: [
              { text: "Sequence ID: ", style: "titleBold" },
              { text: sharedData.sequenceId, fontSize: 11 },
            ],
          },
          {
            text: [
              { text: "Date of Report: ", style: "titleBold" },
              { text: sharedData.dateOfReport, fontSize: 11 },
            ],
          },
          // ===============================================================
          {
            text: [{ text: "Introduction: ", style: "titleBold" }],
          },
          {
            text: [
              {
                text: "This report analyzes the patient's HLA genotype to predict potential adverse drug reactions (ADRs) based on known associations between specific HLA alleles and drug hypersensitivity. The analysis utilizes sequencing reads targeting the HLA-A, HLA-B, and HLA-C genes.",
                style: "normalText",
              },
            ],
          },
          // ===============================================================
          {text:" "},
          {
            text: [{ text: "Methods: ", style: "titleBold" }],
          },
          {
            ol: [
              {
                text: [
                  { text: "HLA Genotyping: ", style: "titleBold" },
                  {
                    text: "Sequencing reads targeting HLA-A, HLA-B, and HLA-C genes were analyzed to determine the patient's HLA genotype.",
                  },
                ],
                style: "normalText",
                margin: [20, 0, 0, 0],
              },
              {
                text: [
                  { text: "ADR Prediction: ", style: "titleBold" },
                  {
                    text: "The identified HLA alleles were cross-referenced with the HLA-ADR database from AlleleFrequencies.net (https://www.allelefrequencies.net/). This database compiles information on HLA allele associations with various ADRs.",
                  },
                ],
                style: "normalText",
                margin: [20, 0, 0, 0],
              },
              {
                text: [
                  { text: "Association Metric: ", style: "titleBold" },
                  {
                    text: "The association metric used consists of an OR (odds ratio) which determines the odds of having an ADR while taking a drug when compared to a normal person and a CI (confidence interval) which are the range on what the true values may possibly be.",
                  },
                ],
                style: "normalText",
                margin: [20, 0, 0, 0],
              },
            ],
          },
          // ===============================================================
          {
            text: [{ text: "Programs: ", style: "titleBold" }],
          },
          {
            ol: [
              {
                text: [
                  { text: "HLA*LA: ", style: "titleBold" },
                  {
                    text: "HLA*LA is a tool designed for high-resolution genotyping of human leukocyte antigen (HLA) genes using next-generation sequencing (NGS) data.",
                  },
                ],
                style: "normalText",
                margin: [20, 0, 0, 0],
              }
            ],
          },
          {text:" "},
          ,
          {
            text: [{ text: "HLA Genotype: ", style: "titleBold" }],
          },
          {
            ul: sharedData.genotypeList,
            style: "normalText",
            margin: [20, 0, 0, 0],
          },
          {text: " "},
          {
            text: [{ text: "Predicted ADRs:", style: "titleBold" }],
          },
          {
            text: [
              {
                text: "Based on the patient's HLA genotype, the following ADRs may be associated with increased risk:",
                style: "normalText",
              },
            ],
          },
          {
            ul: formatADR(),
            margin: [20, 0, 0, 0],
          },
          {text:" "},
          {
            text: [{ text: "Disclaimer: ", style: "titleBold" }],
          },
          {
            text: [
              {
                text: "Discuss these potential ADRs with the patient and their prescribing physician to consider alternative medications if appropriate.",
                style: "normalText",
              },
            ],
          },
          // ===============================================================
          {
            text: [{ text: "References: ", style: "titleBold" }],
          },
          {
            ul: [
              "Ghattaoraya GS, Dundar Y, González-Galarza FF, et al. A web resource for mining HLA associations with adverse drug reactions: HLA-ADR. Database. 2016;2016. https://doi.org/10.1093/database/baw069",
              "Dilthey AT, Mentzer AJ, Carapito R, et al. HLA*LA—HLA typing from linearly projected graph alignments. Bioinformatics. 2019;35(21):4394-4396. https://doi.org/10.1093/bioinformatics/btz235",
            ],
            style: "normalText",
            margin: [20, 0, 0, 0],
          },
          // ===============================================================
          {
            text: [{ text: "Automatically generated by: ", style: "titleBold" }],
          },
          {
            text: [
              { text: "Sira Tongsima with the help of many others.", style: "normalText" },
            ],
          },
        ],
        styles: {
          header: {
            fontSize: 13,
            bold: true,
          },
          titleBold: {
            bold: true,
            fontSize: 12,
          },
          normalText: {
            fontSize: 12,
          },
        },
      };
      pdfMake.fonts = {
        // download default Roboto font from cdnjs.com
        Roboto: {
          normal:
            "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.66/fonts/Roboto/Roboto-Regular.ttf",
          bold: "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.66/fonts/Roboto/Roboto-Medium.ttf",
          italics:
            "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.66/fonts/Roboto/Roboto-Italic.ttf",
          bolditalics:
            "https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.66/fonts/Roboto/Roboto-MediumItalic.ttf",
        },
      };
      // await pdfMake.createPdf(docDefinition).open({}, window);
      await pdfMake.createPdf(docDefinition).download('ADR-Prediction.pdf');
    }
  
  pdfRender().then(() => {router.push("dashboard")});

</script>