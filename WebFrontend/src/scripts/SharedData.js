

function getDate()
{
  const today = new Date();
  const day = String(today.getDate()).padStart(2, '0');
  const month = String(today.getMonth() + 1).padStart(2, '0'); // Months are zero-based
  const year = today.getFullYear();

  const formattedDate = `${day}/${month}/${year}`;
  return formattedDate;
}

export const sharedData = {
  sequenceLabel: "Neymar JR",
  sequenceId: "12345678 2345678 678",
  dateOfReport: getDate(),
  genotypeList: [
    "Null", "Null", "Null", "Null", "Null", "Null",
  ],
  drugList: [
    {allele: "HLA-B*56:01", drug: "oxcarbazepine", odds: "6.16 (0.24 - 156.27)", adr: "maculopapular exanthema (MPE)", link: "http://www.ncbi.nlm.nih.gov/pubmed/23829937"}
  ],
};