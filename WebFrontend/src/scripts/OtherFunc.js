/* eslint-disable*/
function ScrollTo(sectionId) {
      const section = document.getElementById(sectionId);
      if (section) {
        section.scrollIntoView({ behavior: 'smooth' });
      }
}

function ValidName(name)
{
  return name.length > 0;
}

function ValidEmail(email)
{
  const patt = /^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$/;

  return patt.test(email);
}

function ValidPassword(password)
{
  return password.length >= 6;
}

function ValidLabel(label) {
  const patt = /^[^\x00-\x1F\x7F]+$/;
  return patt.test(label) && label.length <= 100;
}

function ErrorHandler(err)
{
  window.sharedData.latestError = err;
  console.error(err);
}

function SilentErrorHandler(err)
{
  console.log("A silent error has occured", err);
}


export{
    ScrollTo,
    ValidName,
    ValidEmail,
    ValidLabel,
    ValidPassword,
    ErrorHandler,
    SilentErrorHandler,
}