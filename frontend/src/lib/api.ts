import axios from 'axios';

// Use empty string for relative path in production/Codespaces, or explicit URL for local dev
const ROOT_API_URL = (import.meta.env.VITE_API_BASE_URL ?? '').replace(
  /\/$/,
  '',
);

const api = axios.create({
  baseURL: `${ROOT_API_URL}/api/v1`,
  timeout: 15000,
});

export default api;
