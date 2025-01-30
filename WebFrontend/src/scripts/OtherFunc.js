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

async function ValidEmail(email)
{
  const patt = /^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$/;

  return patt.test(email);
}

function ValidPassword(password)
{
  return password.length >= 6;
}

function ValidLabel(label)
{
  const patt = /^[a-zA-Z0-9_]+$/;
  return patt.test(label);
}

export{
    ScrollTo,
    ValidName,
    ValidEmail,
    ValidLabel,
    ValidPassword,
}