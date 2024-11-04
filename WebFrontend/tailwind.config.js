/** @type {import('tailwindcss').Config} */
module.exports = {
  content: ['./public/**/*.html', './src/**/*.{vue,js,ts,jsx,tsx}'],
  theme: {
    extend: {
      spacing: {
        '112': '28rem',
        '128': '32rem',
        '144': '36rem',
        '200': '50rem',
        '256': '64rem',
        '512': '128rem'
      },

    },
  },
  daisyui: {
     themes: ["corporate", "sunset"],
  },
  plugins: [
    require('daisyui'),
  ],
}

