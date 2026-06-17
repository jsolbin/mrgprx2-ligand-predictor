import nextVitals from "eslint-config-next/core-web-vitals.js";
import prettierConfig from "eslint-config-prettier";

const config = [
  ...nextVitals,
  prettierConfig,
  {
    ignores: [".next/**", "node_modules/**", "out/**", "coverage/**"],
  },
];

export default config;
