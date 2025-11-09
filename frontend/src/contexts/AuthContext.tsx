import React, { createContext, useContext, useState, useEffect } from 'react';
import api from '../lib/api';

interface User {
  id: number;
  email: string;
  full_name: string;
}

interface AuthContextType {
  user: User | null;
  token: string | null;
  login: (email: string, password: string) => Promise<void>;
  register: (email: string, password: string, fullName: string) => Promise<void>;
  logout: () => void;
  isAuthenticated: boolean;
  isLoading: boolean;
}

const AuthContext = createContext<AuthContextType | undefined>(undefined);

export const AuthProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [user, setUser] = useState<User | null>(null);
  const [token, setToken] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(true);

  useEffect(() => {
    // 从 localStorage 恢复登录状态
    const savedToken = localStorage.getItem('auth_token');
    const savedUser = localStorage.getItem('auth_user');
    
    if (savedToken && savedUser) {
      setToken(savedToken);
      setUser(JSON.parse(savedUser));
      api.defaults.headers.common['Authorization'] = `Bearer ${savedToken}`;
    }
    setIsLoading(false);
  }, []);

  const login = async (email: string, password: string) => {
    try {
      console.log('Attempting login with:', email);
      
      const params = new URLSearchParams();
      params.append('username', email);
      params.append('password', password);

      console.log('Sending request to:', `${api.defaults.baseURL}/auth/login/access-token`);
      
      const response = await api.post('/auth/login/access-token', params, {
        headers: {
          'Content-Type': 'application/x-www-form-urlencoded',
        },
      });

      console.log('Login response:', response.data);
      
      const { access_token } = response.data;
      setToken(access_token);
      
      // 保存 token 和用户信息
      localStorage.setItem('auth_token', access_token);
      api.defaults.headers.common['Authorization'] = `Bearer ${access_token}`;

      // 获取用户信息（简化版本，使用 email）
      const userData: User = {
        id: 1,
        email,
        full_name: email.split('@')[0],
      };
      setUser(userData);
      localStorage.setItem('auth_user', JSON.stringify(userData));
      
      console.log('Login successful!');
    } catch (error: any) {
      console.error('Login error:', error);
      console.error('Error response:', error.response?.data);
      throw error;
    }
  };

  const register = async (email: string, password: string, fullName: string) => {
    try {
      console.log('Attempting registration with:', email);
      
      const response = await api.post('/auth/register', {
        email,
        password,
        full_name: fullName,
      });

      console.log('Registration response:', response.data);
      const userData: User = response.data;
      
      // 注册后自动登录
      console.log('Auto-login after registration...');
      await login(email, password);
    } catch (error: any) {
      console.error('Registration error:', error);
      console.error('Error response:', error.response?.data);
      throw error;
    }
  };

  const logout = () => {
    setUser(null);
    setToken(null);
    localStorage.removeItem('auth_token');
    localStorage.removeItem('auth_user');
    delete api.defaults.headers.common['Authorization'];
  };

  return (
    <AuthContext.Provider
      value={{
        user,
        token,
        login,
        register,
        logout,
        isAuthenticated: !!token,
        isLoading,
      }}
    >
      {children}
    </AuthContext.Provider>
  );
};

export const useAuth = () => {
  const context = useContext(AuthContext);
  if (context === undefined) {
    throw new Error('useAuth must be used within an AuthProvider');
  }
  return context;
};
